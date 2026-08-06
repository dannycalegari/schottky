#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "figure_export.h"

/*****************************************************************************
 * figure_export -- see figure_export.h for what this is and how to use it.
 *
 * Layout of this file:
 *   1. adler32 / crc32
 *   2. DEFLATE (fixed Huffman + LZ77), and the zlib wrapper
 *   3. PNG
 *   4. a 5x7 bitmap font, used only when rasterising text
 *   5. rasterising the vector overlays (PNG, and EPS/PDF with vector_overlays off)
 *   6. EPS
 *   7. PDF
 *   8. the public entry points
 * ***************************************************************************/

namespace figexp {

/*=========================================================== 1. checksums ==*/

unsigned long adler32(const unsigned char* d, size_t n) {
  unsigned long a = 1, b = 0;
  for (size_t i = 0; i < n; ++i) {
    a = (a + d[i]) % 65521;
    b = (b + a) % 65521;
  }
  return (b << 16) | a;
}

static unsigned long crc_table[256];
static bool crc_table_built = false;

static void build_crc_table() {
  for (unsigned long n = 0; n < 256; ++n) {
    unsigned long c = n;
    for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xedb88320UL ^ (c >> 1)) : (c >> 1);
    crc_table[n] = c;
  }
  crc_table_built = true;
}

unsigned long crc32_of(const unsigned char* d, size_t n) {
  if (!crc_table_built) build_crc_table();
  unsigned long c = 0xffffffffUL;
  for (size_t i = 0; i < n; ++i) c = crc_table[(c ^ d[i]) & 0xff] ^ (c >> 8);
  return c ^ 0xffffffffUL;
}

/*============================================================= 2. DEFLATE ==*/

/* A fixed-Huffman DEFLATE compressor.  Fixed Huffman means we do not have to
 * emit a code table, which costs a few percent of ratio and saves a lot of
 * code; on the large flat-color regions typical of these figures the LZ77
 * matcher does nearly all the work anyway.  The point of implementing this at
 * all is to avoid linking zlib -- see the header. */

namespace {

struct BitWriter {
  std::vector<unsigned char>* out;
  unsigned long acc;
  int nbits;
  BitWriter(std::vector<unsigned char>* o) { out = o; acc = 0; nbits = 0; }
  //DEFLATE packs bit fields least-significant-bit first
  void bits(unsigned long v, int n) {
    v &= (n >= 32 ? 0xffffffffUL : ((1UL << n) - 1));
    acc |= (v << nbits);
    nbits += n;
    while (nbits >= 8) {
      out->push_back((unsigned char)(acc & 0xff));
      acc >>= 8;
      nbits -= 8;
    }
  }
  //...but Huffman codes are packed most-significant-bit first
  void huff(unsigned long code, int n) {
    for (int i = n - 1; i >= 0; --i) bits((code >> i) & 1, 1);
  }
  void flush() {
    if (nbits > 0) {
      out->push_back((unsigned char)(acc & 0xff));
      acc = 0;
      nbits = 0;
    }
  }
};

//RFC 1951 section 3.2.5
const int LEN_BASE[29]  = {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,
                           59,67,83,99,115,131,163,195,227,258};
const int LEN_EXTRA[29] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                           5,5,5,5,0};
const int DIST_BASE[30] = {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,
                           513,769,1025,1537,2049,3073,4097,6145,8193,12289,
                           16385,24577};
const int DIST_EXTRA[30]= {0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
                           11,11,12,12,13,13};

//RFC 1951 section 3.2.6: the fixed literal/length code
void fixed_lit(int sym, unsigned long* code, int* len) {
  if (sym < 144)      { *code = 0x30UL  + sym;         *len = 8; }
  else if (sym < 256) { *code = 0x190UL + (sym - 144); *len = 9; }
  else if (sym < 280) { *code = (unsigned long)(sym - 256); *len = 7; }
  else                { *code = 0xC0UL  + (sym - 280); *len = 8; }
}

inline void emit_literal(BitWriter& bw, int byte) {
  unsigned long c; int l;
  fixed_lit(byte, &c, &l);
  bw.huff(c, l);
}

void emit_match(BitWriter& bw, int length, int dist) {
  int li = 28;
  while (li > 0 && LEN_BASE[li] > length) --li;
  unsigned long c; int l;
  fixed_lit(257 + li, &c, &l);
  bw.huff(c, l);
  if (LEN_EXTRA[li]) bw.bits((unsigned long)(length - LEN_BASE[li]), LEN_EXTRA[li]);

  int di = 29;
  while (di > 0 && DIST_BASE[di] > dist) --di;
  bw.huff((unsigned long)di, 5);                //distance codes are 5 bits fixed
  if (DIST_EXTRA[di]) bw.bits((unsigned long)(dist - DIST_BASE[di]), DIST_EXTRA[di]);
}

const int WINDOW   = 32768;
const int MIN_MATCH = 3;
const int MAX_MATCH = 258;
const int HASH_BITS = 15;
const int HASH_SIZE = 1 << HASH_BITS;
const int MAX_CHAIN = 128;      //how hard to look for a longer match

inline unsigned int hash3(const unsigned char* p) {
  unsigned long v = ((unsigned long)p[0] << 16) | ((unsigned long)p[1] << 8) | p[2];
  return (unsigned int)((v * 2654435761UL) >> (32 - HASH_BITS)) & (HASH_SIZE - 1);
}

} //anonymous namespace

namespace {

/* A DEFLATE stream of stored (uncompressed) blocks.  Fixed Huffman costs 8 or 9
 * bits per literal, so on incompressible data it EXPANDS by a few percent; a
 * stored block never expands by more than 5 bytes per 65535.  deflate_raw
 * compresses both ways and keeps the smaller, which bounds the output. */
void deflate_stored(const unsigned char* d, size_t n, std::vector<unsigned char>& out) {
  size_t i = 0;
  do {
    size_t chunk = n - i;
    if (chunk > 65535) chunk = 65535;
    bool last = (i + chunk >= n);
    out.push_back(last ? 1 : 0);                //BFINAL in bit 0, BTYPE 00
    out.push_back((unsigned char)( chunk        & 0xff));
    out.push_back((unsigned char)((chunk >> 8)  & 0xff));
    out.push_back((unsigned char)((~chunk)      & 0xff));
    out.push_back((unsigned char)((~chunk >> 8) & 0xff));
    for (size_t k = 0; k < chunk; ++k) out.push_back(d[i + k]);
    i += chunk;
  } while (i < n);
}

} //anonymous namespace

void deflate_raw(const unsigned char* d, size_t n, std::vector<unsigned char>& out) {
  std::vector<unsigned char> huff_out;
  BitWriter bw(&huff_out);
  bw.bits(1, 1);                                //BFINAL = 1, one block
  bw.bits(1, 2);                                //BTYPE  = 01, fixed Huffman

  std::vector<int> head(HASH_SIZE, -1);
  std::vector<int> prev(n > 0 ? n : 1, -1);

  size_t i = 0;
  while (i < n) {
    int best_len = 0, best_dist = 0;
    if (i + MIN_MATCH <= n) {
      unsigned int h = hash3(d + i);
      int j = head[h];
      int chain = 0;
      int maxl = (int)(n - i);
      if (maxl > MAX_MATCH) maxl = MAX_MATCH;
      while (j >= 0 && (int)(i - j) <= WINDOW && chain < MAX_CHAIN) {
        int l = 0;
        while (l < maxl && d[j + l] == d[i + l]) ++l;
        if (l > best_len) {
          best_len = l;
          best_dist = (int)(i - j);
          if (l >= maxl) break;
        }
        j = prev[j];
        ++chain;
      }
      prev[i] = head[h];
      head[h] = (int)i;
    }

    if (best_len >= MIN_MATCH) {
      emit_match(bw, best_len, best_dist);
      //index the interior of the match too, so later matches can find it
      for (size_t k = i + 1; k < i + (size_t)best_len; ++k) {
        if (k + MIN_MATCH > n) break;
        unsigned int h2 = hash3(d + k);
        prev[k] = head[h2];
        head[h2] = (int)k;
      }
      i += (size_t)best_len;
    } else {
      emit_literal(bw, d[i]);
      ++i;
    }
  }

  unsigned long c; int l;
  fixed_lit(256, &c, &l);                       //end of block
  bw.huff(c, l);
  bw.flush();

  //keep whichever encoding is smaller (see deflate_stored)
  std::vector<unsigned char> stored;
  deflate_stored(d, n, stored);
  const std::vector<unsigned char>& best = (stored.size() < huff_out.size() ? stored : huff_out);
  out.insert(out.end(), best.begin(), best.end());
}

void zlib_wrap(const unsigned char* d, size_t n, std::vector<unsigned char>& out) {
  out.push_back(0x78);                          //CM = 8 (deflate), CINFO = 7 (32K)
  out.push_back(0x9c);                          //FLEVEL = 2, FCHECK making it %31==0
  deflate_raw(d, n, out);
  unsigned long a = adler32(d, n);
  out.push_back((unsigned char)((a >> 24) & 0xff));
  out.push_back((unsigned char)((a >> 16) & 0xff));
  out.push_back((unsigned char)((a >>  8) & 0xff));
  out.push_back((unsigned char)( a        & 0xff));
}

/*================================================================== 3. PNG ==*/

namespace {

void be32(std::vector<unsigned char>& v, unsigned long x) {
  v.push_back((unsigned char)((x >> 24) & 0xff));
  v.push_back((unsigned char)((x >> 16) & 0xff));
  v.push_back((unsigned char)((x >>  8) & 0xff));
  v.push_back((unsigned char)( x        & 0xff));
}

void png_chunk(std::vector<unsigned char>& f, const char* tag,
               const std::vector<unsigned char>& data) {
  be32(f, (unsigned long)data.size());
  std::vector<unsigned char> td;
  td.reserve(4 + data.size());
  for (int i = 0; i < 4; ++i) td.push_back((unsigned char)tag[i]);
  td.insert(td.end(), data.begin(), data.end());
  f.insert(f.end(), td.begin(), td.end());
  be32(f, crc32_of(td.empty() ? (const unsigned char*)"" : &td[0], td.size()));
}

//PNG scanline filters.  Choosing per row by least sum of absolute values is the
//standard cheap heuristic and is worth a lot on figures with flat regions.
void filter_row(const unsigned char* cur, const unsigned char* up, int nbytes,
                int bpp, std::vector<unsigned char>& out,
                std::vector<unsigned char>* cand) {
  long best_score = -1;
  int best = 0;
  for (int f = 0; f < 3; ++f) {
    cand[f].resize(nbytes);
    long score = 0;
    for (int k = 0; k < nbytes; ++k) {
      int a = (k >= bpp) ? cur[k - bpp] : 0;
      int b = up ? up[k] : 0;
      int v;
      if (f == 0)      v = cur[k];
      else if (f == 1) v = cur[k] - a;
      else             v = cur[k] - b;
      cand[f][k] = (unsigned char)(v & 0xff);
      int s = (signed char)cand[f][k];
      score += (s < 0 ? -s : s);
    }
    if (best_score < 0 || score < best_score) { best_score = score; best = f; }
  }
  out.push_back((unsigned char)best);
  out.insert(out.end(), cand[best].begin(), cand[best].end());
}

} //anonymous namespace

/* The workhorse.  `trns` asks for a tRNS chunk naming one fully transparent
 * color, which is how a truecolor PNG expresses transparency without carrying
 * a whole alpha channel. */
static bool write_png_impl(const std::string& path, int w, int h,
                           const std::vector<unsigned char>& rgb_bottom_up,
                           bool trns, unsigned char tr, unsigned char tg,
                           unsigned char tb, std::string* err) {
  if (w <= 0 || h <= 0 || (int)rgb_bottom_up.size() < 3*w*h) {
    if (err) *err = "write_png: bad size or short buffer";
    return false;
  }
  //PNG rows run top to bottom; our buffer runs bottom to top
  std::vector<unsigned char> raw;
  raw.reserve((size_t)h * (3*(size_t)w + 1));
  std::vector<unsigned char> cand[3];               //scratch, reused per row
  for (int y = h - 1; y >= 0; --y) {
    const unsigned char* cur = &rgb_bottom_up[3*(size_t)y*w];
    //"the row above" means the previously written one, and we write top-down
    const unsigned char* up  = (y == h - 1) ? 0 : &rgb_bottom_up[3*(size_t)(y+1)*w];
    filter_row(cur, up, 3*w, 3, raw, cand);
  }
  std::vector<unsigned char> idat;
  zlib_wrap(raw.empty() ? (const unsigned char*)"" : &raw[0], raw.size(), idat);

  std::vector<unsigned char> f;
  const unsigned char sig[8] = {0x89,'P','N','G',0x0d,0x0a,0x1a,0x0a};
  f.insert(f.end(), sig, sig + 8);

  std::vector<unsigned char> ihdr;
  be32(ihdr, (unsigned long)w);
  be32(ihdr, (unsigned long)h);
  ihdr.push_back(8);      //bit depth
  ihdr.push_back(2);      //color type 2 = truecolor RGB
  ihdr.push_back(0);      //deflate
  ihdr.push_back(0);      //adaptive filtering
  ihdr.push_back(0);      //no interlace
  png_chunk(f, "IHDR", ihdr);
  if (trns) {
    //tRNS for color type 2 is one RGB triple, each channel a 16-bit big-endian
    //sample; every pixel of exactly that color becomes fully transparent
    std::vector<unsigned char> t;
    t.push_back(0); t.push_back(tr);
    t.push_back(0); t.push_back(tg);
    t.push_back(0); t.push_back(tb);
    png_chunk(f, "tRNS", t);
  }
  png_chunk(f, "IDAT", idat);
  png_chunk(f, "IEND", std::vector<unsigned char>());

  std::ofstream fp(path.c_str(), std::ios::binary);
  if (!fp) { if (err) *err = "cannot open " + path + " for writing"; return false; }
  fp.write((const char*)&f[0], (std::streamsize)f.size());
  fp.flush();
  if (!fp) { if (err) *err = "write failed on " + path; return false; }
  return true;
}

bool write_png(const std::string& path, int w, int h,
               const std::vector<unsigned char>& rgb_bottom_up, std::string* err) {
  return write_png_impl(path, w, h, rgb_bottom_up, false, 0, 0, 0, err);
}

/*=========================================================== 4. 5x7 font ===*/

/* Used only when rasterising text into a PNG.  EPS and PDF use Helvetica,
 * which is built into both formats and is what a paper figure wants.  Five
 * columns per glyph, one byte each, bit 0 = top row. */

namespace {

const unsigned char FONT5x7[95][5] = {
  {0x00,0x00,0x00,0x00,0x00}, /*   */ {0x00,0x00,0x5F,0x00,0x00}, /* ! */
  {0x00,0x07,0x00,0x07,0x00}, /* " */ {0x14,0x7F,0x14,0x7F,0x14}, /* # */
  {0x24,0x2A,0x7F,0x2A,0x12}, /* $ */ {0x23,0x13,0x08,0x64,0x62}, /* % */
  {0x36,0x49,0x55,0x22,0x50}, /* & */ {0x00,0x05,0x03,0x00,0x00}, /* ' */
  {0x00,0x1C,0x22,0x41,0x00}, /* ( */ {0x00,0x41,0x22,0x1C,0x00}, /* ) */
  {0x14,0x08,0x3E,0x08,0x14}, /* * */ {0x08,0x08,0x3E,0x08,0x08}, /* + */
  {0x00,0x50,0x30,0x00,0x00}, /* , */ {0x08,0x08,0x08,0x08,0x08}, /* - */
  {0x00,0x60,0x60,0x00,0x00}, /* . */ {0x20,0x10,0x08,0x04,0x02}, /* / */
  {0x3E,0x51,0x49,0x45,0x3E}, /* 0 */ {0x00,0x42,0x7F,0x40,0x00}, /* 1 */
  {0x42,0x61,0x51,0x49,0x46}, /* 2 */ {0x21,0x41,0x45,0x4B,0x31}, /* 3 */
  {0x18,0x14,0x12,0x7F,0x10}, /* 4 */ {0x27,0x45,0x45,0x45,0x39}, /* 5 */
  {0x3C,0x4A,0x49,0x49,0x30}, /* 6 */ {0x01,0x71,0x09,0x05,0x03}, /* 7 */
  {0x36,0x49,0x49,0x49,0x36}, /* 8 */ {0x06,0x49,0x49,0x29,0x1E}, /* 9 */
  {0x00,0x36,0x36,0x00,0x00}, /* : */ {0x00,0x56,0x36,0x00,0x00}, /* ; */
  {0x08,0x14,0x22,0x41,0x00}, /* < */ {0x14,0x14,0x14,0x14,0x14}, /* = */
  {0x00,0x41,0x22,0x14,0x08}, /* > */ {0x02,0x01,0x51,0x09,0x06}, /* ? */
  {0x32,0x49,0x79,0x41,0x3E}, /* @ */ {0x7E,0x11,0x11,0x11,0x7E}, /* A */
  {0x7F,0x49,0x49,0x49,0x36}, /* B */ {0x3E,0x41,0x41,0x41,0x22}, /* C */
  {0x7F,0x41,0x41,0x22,0x1C}, /* D */ {0x7F,0x49,0x49,0x49,0x41}, /* E */
  {0x7F,0x09,0x09,0x09,0x01}, /* F */ {0x3E,0x41,0x49,0x49,0x7A}, /* G */
  {0x7F,0x08,0x08,0x08,0x7F}, /* H */ {0x00,0x41,0x7F,0x41,0x00}, /* I */
  {0x20,0x40,0x41,0x3F,0x01}, /* J */ {0x7F,0x08,0x14,0x22,0x41}, /* K */
  {0x7F,0x40,0x40,0x40,0x40}, /* L */ {0x7F,0x02,0x0C,0x02,0x7F}, /* M */
  {0x7F,0x04,0x08,0x10,0x7F}, /* N */ {0x3E,0x41,0x41,0x41,0x3E}, /* O */
  {0x7F,0x09,0x09,0x09,0x06}, /* P */ {0x3E,0x41,0x51,0x21,0x5E}, /* Q */
  {0x7F,0x09,0x19,0x29,0x46}, /* R */ {0x46,0x49,0x49,0x49,0x31}, /* S */
  {0x01,0x01,0x7F,0x01,0x01}, /* T */ {0x3F,0x40,0x40,0x40,0x3F}, /* U */
  {0x1F,0x20,0x40,0x20,0x1F}, /* V */ {0x7F,0x20,0x18,0x20,0x7F}, /* W */
  {0x63,0x14,0x08,0x14,0x63}, /* X */ {0x03,0x04,0x78,0x04,0x03}, /* Y */
  {0x61,0x51,0x49,0x45,0x43}, /* Z */ {0x00,0x7F,0x41,0x41,0x00}, /* [ */
  {0x02,0x04,0x08,0x10,0x20}, /* \ */ {0x00,0x41,0x41,0x7F,0x00}, /* ] */
  {0x04,0x02,0x01,0x02,0x04}, /* ^ */ {0x40,0x40,0x40,0x40,0x40}, /* _ */
  {0x00,0x01,0x02,0x04,0x00}, /* ` */ {0x20,0x54,0x54,0x54,0x78}, /* a */
  {0x7F,0x48,0x44,0x44,0x38}, /* b */ {0x38,0x44,0x44,0x44,0x20}, /* c */
  {0x38,0x44,0x44,0x48,0x7F}, /* d */ {0x38,0x54,0x54,0x54,0x18}, /* e */
  {0x08,0x7E,0x09,0x01,0x02}, /* f */ {0x0C,0x52,0x52,0x52,0x3E}, /* g */
  {0x7F,0x08,0x04,0x04,0x78}, /* h */ {0x00,0x44,0x7D,0x40,0x00}, /* i */
  {0x20,0x40,0x44,0x3D,0x00}, /* j */ {0x7F,0x10,0x28,0x44,0x00}, /* k */
  {0x00,0x41,0x7F,0x40,0x00}, /* l */ {0x7C,0x04,0x18,0x04,0x78}, /* m */
  {0x7C,0x08,0x04,0x04,0x78}, /* n */ {0x38,0x44,0x44,0x44,0x38}, /* o */
  {0x7C,0x14,0x14,0x14,0x08}, /* p */ {0x08,0x14,0x14,0x18,0x7C}, /* q */
  {0x7C,0x08,0x04,0x04,0x08}, /* r */ {0x48,0x54,0x54,0x54,0x20}, /* s */
  {0x04,0x3F,0x44,0x40,0x20}, /* t */ {0x3C,0x40,0x40,0x20,0x1C}, /* u */
  {0x1C,0x20,0x40,0x20,0x1C}, /* v */ {0x3C,0x40,0x30,0x40,0x3C}, /* w */
  {0x44,0x28,0x10,0x28,0x44}, /* x */ {0x0C,0x50,0x50,0x50,0x3C}, /* y */
  {0x44,0x64,0x54,0x4C,0x44}, /* z */ {0x00,0x08,0x36,0x41,0x00}, /* { */
  {0x00,0x00,0x7F,0x00,0x00}, /* | */ {0x00,0x41,0x36,0x08,0x00}, /* } */
  {0x08,0x04,0x08,0x10,0x08}  /* ~ */
};

//advance per character, in glyph cells (5 columns + 1 of spacing)
const int FONT_ADV = 6;
const int FONT_W = 5, FONT_H = 7;

//Helvetica's cap height as a fraction of the em; see Canvas::text
const double FONT_CAP_HEIGHT_EM = 0.717;

/* Helvetica character widths, in 1/1000 em, for ASCII 32..126 (from the
 * standard AFM).  PDF, unlike PostScript, has no stringwidth operator, so
 * centring and right-aligning text in the PDF back end needs the metrics
 * here; otherwise anchored labels drift by several percent of their length. */
const short HELVETICA_W[95] = {
  278,278,355,556,556,889,667,191,333,333,389,584,278,333,278,278,  /*  sp../  */
  556,556,556,556,556,556,556,556,556,556,                          /*  0..9   */
  278,278,584,584,584,556,1015,                                     /*  :..@   */
  667,667,722,722,667,611,778,722,278,500,667,556,833,722,778,      /*  A..O   */
  667,778,722,667,611,722,667,944,667,667,611,                      /*  P..Z   */
  278,278,278,469,556,333,                                          /*  [..`   */
  556,556,500,556,556,278,556,556,222,222,500,222,833,556,556,      /*  a..o   */
  556,556,333,500,278,556,500,722,500,500,500,                      /*  p..z   */
  334,260,334,584                                                   /*  {..~   */
};

double helvetica_width(const std::string& s, double size_pt) {
  double w = 0.0;
  for (size_t i = 0; i < s.size(); ++i) {
    int c = (unsigned char)s[i];
    if (c < 32 || c > 126) c = '?';
    w += HELVETICA_W[c - 32];
  }
  return w*size_pt/1000.0;
}

} //anonymous namespace

/*=================================================== 5. rasterising paths ==*/

namespace {

struct Xform {
  //mathematical coordinates -> device coordinates (y up in both)
  double x0, y0, sx, sy;
  double to_x(double x) const { return (x - x0) * sx; }
  double to_y(double y) const { return (y - y0) * sy; }
};

struct Canvas {
  Raster* R;
  Xform T;
  double px_per_pt;

  /* Clip a device-space segment to a margin around the canvas (Liang-Barsky).
   * Two reasons this matters rather than relying on the per-pixel bounds check
   * in blend(): a segment a million pixels long would step through it half a
   * pixel at a time, and a coordinate that came out as inf or nan -- easy
   * enough if a caller hands us an unset point -- would loop forever.  Returns
   * false if the segment misses the canvas entirely. */
  bool clip_segment(double* ax, double* ay, double* bx, double* by, double margin) const {
    //reject non-finite input outright; there is nothing sensible to draw
    if (!(*ax == *ax) || !(*ay == *ay) || !(*bx == *bx) || !(*by == *by)) return false;
    const double LIM = 1e12;
    if (*ax < -LIM || *ax > LIM || *ay < -LIM || *ay > LIM) return false;
    if (*bx < -LIM || *bx > LIM || *by < -LIM || *by > LIM) return false;

    double xlo = -margin, xhi = R->w + margin;
    double ylo = -margin, yhi = R->h + margin;
    double dx = *bx - *ax, dy = *by - *ay;
    double t0 = 0.0, t1 = 1.0;
    const double p[4] = {-dx, dx, -dy, dy};
    const double q[4] = {*ax - xlo, xhi - *ax, *ay - ylo, yhi - *ay};
    for (int i = 0; i < 4; ++i) {
      if (p[i] == 0.0) {
        if (q[i] < 0.0) return false;             //parallel and outside
      } else {
        double r = q[i]/p[i];
        if (p[i] < 0.0) { if (r > t1) return false; if (r > t0) t0 = r; }
        else            { if (r < t0) return false; if (r < t1) t1 = r; }
      }
    }
    double nax = *ax + t0*dx, nay = *ay + t0*dy;
    double nbx = *ax + t1*dx, nby = *ay + t1*dy;
    *ax = nax; *ay = nay; *bx = nbx; *by = nby;
    return true;
  }

  void blend(int i, int j, const Style& st, double cov) {
    if (i < 0 || j < 0 || i >= R->w || j >= R->h) return;
    if (cov <= 0) return;
    if (cov > 1) cov = 1;
    unsigned char r, g, b;
    R->get_pixel(i, j, &r, &g, &b);
    double nr = r*(1-cov) + 255.0*st.r*cov;
    double ng = g*(1-cov) + 255.0*st.g*cov;
    double nb = b*(1-cov) + 255.0*st.b*cov;
    R->set_pixel(i, j, (unsigned char)(nr + 0.5),
                       (unsigned char)(ng + 0.5),
                       (unsigned char)(nb + 0.5));
  }

  //a filled disk in device pixels, antialiased at the rim
  void disk_px(double cx, double cy, double rad, const Style& st) {
    if (!(cx == cx) || !(cy == cy) || !(rad == rad)) return;      //nan
    if (rad < 0.35) rad = 0.35;
    int i0 = (int)std::floor(cx - rad - 1), i1 = (int)std::ceil(cx + rad + 1);
    int j0 = (int)std::floor(cy - rad - 1), j1 = (int)std::ceil(cy + rad + 1);
    //clamp to the canvas: blend() would reject the pixels anyway, but a disk of
    //radius 1e9 must not cost 1e18 iterations to reject
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 > R->w - 1) i1 = R->w - 1;
    if (j1 > R->h - 1) j1 = R->h - 1;
    for (int j = j0; j <= j1; ++j) {
      for (int i = i0; i <= i1; ++i) {
        double dx = (i + 0.5) - cx, dy = (j + 0.5) - cy;
        double d = std::sqrt(dx*dx + dy*dy);
        blend(i, j, st, rad + 0.5 - d);          //1px feather
      }
    }
  }

  void segment_px(double ax, double ay, double bx, double by, double lw,
                  const Style& st) {
    double r = lw/2.0;
    if (!clip_segment(&ax, &ay, &bx, &by, r + 2.0)) return;
    double dx = bx - ax, dy = by - ay;
    double len = std::sqrt(dx*dx + dy*dy);
    if (len < 1e-9) { disk_px(ax, ay, r, st); return; }
    int steps = (int)(len*2.0) + 1;
    for (int k = 0; k <= steps; ++k) {
      double t = double(k)/steps;
      disk_px(ax + dx*t, ay + dy*t, r, st);
    }
  }

  //one segment of a dashed line.  phase carries the arc length already travelled
  //along the whole polyline, so dashes continue across corners as they do in
  //PostScript and PDF.
  void dashed_segment_px(double ax, double ay, double bx, double by, double lw,
                         const Style& st, double* phase) {
    double on = 3*lw, off = 2*lw, period = on + off;
    if (period < 1e-6) { segment_px(ax, ay, bx, by, lw, st); return; }
    if (!clip_segment(&ax, &ay, &bx, &by, lw + 2.0)) return;
    double dx = bx - ax, dy = by - ay;
    double len = std::sqrt(dx*dx + dy*dy);
    if (len < 1e-9) return;
    double t = 0.0;
    const double step = 0.5;                   //half a pixel at a time
    while (t < len) {
      double t2 = t + step;
      if (t2 > len) t2 = len;
      double ph = std::fmod(*phase + t, period);
      if (ph < on) {
        double u = t/len, v = t2/len;
        segment_px(ax + dx*u, ay + dy*u, ax + dx*v, ay + dy*v, lw, st);
      }
      t = t2;
    }
    *phase = std::fmod(*phase + len, period);
  }

  void polyline(const std::vector<Point2d<double> >& p, bool closed, double lw,
                const Style& st) {
    if (p.size() < 2) {
      if (p.size() == 1) disk_px(T.to_x(p[0].x), T.to_y(p[0].y), lw/2.0, st);
      return;
    }
    double phase = 0.0;
    for (size_t k = 0; k + 1 < p.size(); ++k) {
      double ax = T.to_x(p[k].x),   ay = T.to_y(p[k].y);
      double bx = T.to_x(p[k+1].x), by = T.to_y(p[k+1].y);
      if (st.dashed) dashed_segment_px(ax, ay, bx, by, lw, st, &phase);
      else           segment_px(ax, ay, bx, by, lw, st);
    }
    if (closed) {
      double ax = T.to_x(p.back().x), ay = T.to_y(p.back().y);
      double bx = T.to_x(p[0].x),     by = T.to_y(p[0].y);
      if (st.dashed) dashed_segment_px(ax, ay, bx, by, lw, st, &phase);
      else           segment_px(ax, ay, bx, by, lw, st);
    }
  }

  //even-odd scanline fill of a polygon given in mathematical coordinates
  void fill_polygon(const std::vector<Point2d<double> >& p, const Style& st) {
    if (p.size() < 3) return;
    std::vector<double> X(p.size()), Y(p.size());
    double ymin = 1e300, ymax = -1e300;
    for (size_t k = 0; k < p.size(); ++k) {
      X[k] = T.to_x(p[k].x);
      Y[k] = T.to_y(p[k].y);
      if (Y[k] < ymin) ymin = Y[k];
      if (Y[k] > ymax) ymax = Y[k];
    }
    int j0 = (int)std::floor(ymin), j1 = (int)std::ceil(ymax);
    if (j0 < 0) j0 = 0;
    if (j1 > R->h - 1) j1 = R->h - 1;
    std::vector<double> xs;
    for (int j = j0; j <= j1; ++j) {
      double yc = j + 0.5;
      xs.clear();
      for (size_t k = 0; k < p.size(); ++k) {
        size_t k2 = (k + 1) % p.size();
        double ya = Y[k], yb = Y[k2];
        if ((ya <= yc && yb > yc) || (yb <= yc && ya > yc)) {
          double t = (yc - ya)/(yb - ya);
          xs.push_back(X[k] + t*(X[k2] - X[k]));
        }
      }
      if (xs.size() < 2) continue;
      std::sort(xs.begin(), xs.end());
      for (size_t k = 0; k + 1 < xs.size(); k += 2) {
        double xa = xs[k], xb = xs[k+1];
        if (xa < 0.0) xa = 0.0;                       //clamp, so a span reaching
        if (xb > R->w) xb = R->w;                     //to 1e9 is not iterated
        if (xb <= xa) continue;
        int ia = (int)std::floor(xa + 0.5), ib = (int)std::floor(xb + 0.5);
        for (int i = ia; i < ib; ++i) blend(i, j, st, 1.0);
      }
    }
  }

  void text(double mx, double my, const std::string& s, double size_pt, int anchor,
            const Style& st) {
    //Scale so that the 7-row glyph box is the CAP HEIGHT of the nominal point
    //size, not the whole em.  Helvetica -- which the EPS and PDF back ends use
    //at the same nominal size -- has a cap height near 0.72 em, so without this
    //factor the same Figure would come out with visibly larger text in PNG.
    double scale = (FONT_CAP_HEIGHT_EM * size_pt * px_per_pt) / double(FONT_H);
    if (scale < 1.0) scale = 1.0;
    double advance = FONT_ADV * scale;
    double total = advance * (double)s.size();
    double x = T.to_x(mx);
    double y = T.to_y(my);
    if (anchor == 1) x -= total/2.0;
    else if (anchor == 2) x -= total;
    for (size_t c = 0; c < s.size(); ++c) {
      int ch = (unsigned char)s[c];
      if (ch < 32 || ch > 126) ch = '?';
      const unsigned char* g = FONT5x7[ch - 32];
      for (int col = 0; col < FONT_W; ++col) {
        for (int row = 0; row < FONT_H; ++row) {
          if (!((g[col] >> row) & 1)) continue;
          //font rows run downwards, device y runs up
          double gx = x + (col + 0.5)*scale;
          double gy = y + (FONT_H - 1 - row + 0.5)*scale;
          double half = scale/2.0;
          int i0 = (int)std::floor(gx - half), i1 = (int)std::ceil(gx + half);
          int j0 = (int)std::floor(gy - half), j1 = (int)std::ceil(gy + half);
          for (int j = j0; j < j1; ++j)
            for (int i = i0; i < i1; ++i) blend(i, j, st, 1.0);
        }
      }
      x += advance;
    }
  }
};

//expand a Disk into a polygon, so the rasteriser only needs polylines
void circle_pts(double cx, double cy, double r, int nseg,
                std::vector<Point2d<double> >& out) {
  out.clear();
  for (int k = 0; k < nseg; ++k) {
    double th = 2.0*M_PI*k/nseg;
    out.push_back(Point2d<double>(cx + r*std::cos(th), cy + r*std::sin(th)));
  }
}

void rasterise_overlays(const Figure& F, const Options& opt, Canvas& C) {
  double lw_default = opt.line_width_pt * C.px_per_pt;
  for (size_t k = 0; k < F.paths.size(); ++k) {
    const Path& p = F.paths[k];
    double lw = (p.st.line_width > 0 ? p.st.line_width : opt.line_width_pt) * C.px_per_pt;
    if (p.st.filled) C.fill_polygon(p.pts, p.st);
    else             C.polyline(p.pts, p.closed, lw, p.st);
  }
  std::vector<Point2d<double> > cp;
  for (size_t k = 0; k < F.disks.size(); ++k) {
    const Disk& d = F.disks[k];
    circle_pts(d.c.x, d.c.y, d.r, 128, cp);
    double lw = (d.st.line_width > 0 ? d.st.line_width : opt.line_width_pt) * C.px_per_pt;
    if (d.st.filled) C.fill_polygon(cp, d.st);
    else             C.polyline(cp, true, lw, d.st);
  }
  for (size_t k = 0; k < F.dots.size(); ++k) {
    const Dot& d = F.dots[k];
    C.disk_px(C.T.to_x(d.c.x), C.T.to_y(d.c.y), d.r_pt * C.px_per_pt, d.st);
  }
  for (size_t k = 0; k < F.texts.size(); ++k) {
    const Text& t = F.texts[k];
    C.text(t.p.x, t.p.y, t.s, t.size_pt, t.anchor, t.st);
  }
  (void)lw_default;
}

//scale a raster to a target width by nearest neighbour (the source is itself a
//sampling of a continuum, so smoothing it would only blur certain structure)
void resample(const Raster& src, int W, int H, Raster& dst) {
  dst.set_size(W, H);
  for (int j = 0; j < H; ++j) {
    int sj = (int)((double)j*src.h/H);
    if (sj >= src.h) sj = src.h - 1;
    for (int i = 0; i < W; ++i) {
      int si = (int)((double)i*src.w/W);
      if (si >= src.w) si = src.w - 1;
      unsigned char r, g, b;
      src.get_pixel(si, sj, &r, &g, &b);
      dst.set_pixel(i, j, r, g, b);
    }
  }
}

//the frame and axes, added as ordinary overlays so every back end gets them
void add_decorations(Figure& F, const Options& opt) {
  if (opt.draw_axes) {
    Style ax = Style::stroke(0.55, 0.55, 0.55, 0.4);
    if (F.x0 < 0.0 && F.x1 > 0.0) F.add_segment(0.0, F.y0, 0.0, F.y1, ax);
    if (F.y0 < 0.0 && F.y1 > 0.0) F.add_segment(F.x0, 0.0, F.x1, 0.0, ax);
  }
  if (opt.draw_frame) {
    Style fr = Style::stroke(0, 0, 0, 0.5);
    std::vector<Point2d<double> > b;
    b.push_back(Point2d<double>(F.x0, F.y0));
    b.push_back(Point2d<double>(F.x1, F.y0));
    b.push_back(Point2d<double>(F.x1, F.y1));
    b.push_back(Point2d<double>(F.x0, F.y1));
    F.add_path(b, fr, true);
  }
}

} //anonymous namespace

/*================================================================== 6. EPS ==*/

namespace {

std::string ps_escape(const std::string& s) {
  std::string o;
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    if (c == '(' || c == ')' || c == '\\') o += '\\';
    o += c;
  }
  return o;
}

void ps_style(std::ostringstream& o, const Style& st, double lw_default) {
  o << st.r << " " << st.g << " " << st.b << " setrgbcolor\n";
  double lw = (st.line_width > 0 ? st.line_width : lw_default);
  o << lw << " setlinewidth\n";
  if (st.dashed) o << "[" << 3*lw << " " << 2*lw << "] 0 setdash\n";
  else           o << "[] 0 setdash\n";
}

void ps_emit_overlays(std::ostringstream& o, const Figure& F, const Options& opt,
                      const Xform& T) {
  for (size_t k = 0; k < F.paths.size(); ++k) {
    const Path& p = F.paths[k];
    if (p.pts.size() < 2) continue;
    ps_style(o, p.st, opt.line_width_pt);
    o << T.to_x(p.pts[0].x) << " " << T.to_y(p.pts[0].y) << " M\n";
    for (size_t m = 1; m < p.pts.size(); ++m)
      o << T.to_x(p.pts[m].x) << " " << T.to_y(p.pts[m].y) << " L\n";
    if (p.closed) o << "closepath\n";
    o << (p.st.filled ? "fill\n" : "stroke\n");
  }
  for (size_t k = 0; k < F.disks.size(); ++k) {
    const Disk& d = F.disks[k];
    ps_style(o, d.st, opt.line_width_pt);
    o << "newpath " << T.to_x(d.c.x) << " " << T.to_y(d.c.y) << " "
      << d.r*T.sx << " 0 360 arc " << (d.st.filled ? "fill\n" : "stroke\n");
  }
  for (size_t k = 0; k < F.dots.size(); ++k) {
    const Dot& d = F.dots[k];
    ps_style(o, d.st, opt.line_width_pt);
    o << "newpath " << T.to_x(d.c.x) << " " << T.to_y(d.c.y) << " "
      << d.r_pt << " 0 360 arc fill\n";
  }
  for (size_t k = 0; k < F.texts.size(); ++k) {
    const Text& t = F.texts[k];
    o << t.st.r << " " << t.st.g << " " << t.st.b << " setrgbcolor\n";
    o << "/Helvetica findfont " << t.size_pt << " scalefont setfont\n";
    std::string s = ps_escape(t.s);
    o << T.to_x(t.p.x) << " " << T.to_y(t.p.y) << " M\n";
    if (t.anchor == 0)      o << "(" << s << ") show\n";
    else if (t.anchor == 1) o << "(" << s << ") dup stringwidth pop 2 div neg 0 rmoveto show\n";
    else                    o << "(" << s << ") dup stringwidth pop neg 0 rmoveto show\n";
  }
}

bool write_eps(const Figure& F, const Options& opt, const Raster& img,
               double Wpt, double Hpt, const std::string& path, std::string* err) {
  Xform T;
  T.x0 = F.x0; T.y0 = F.y0;
  T.sx = Wpt/(F.x1 - F.x0);
  T.sy = Hpt/(F.y1 - F.y0);

  std::ostringstream o;
  o.setf(std::ios::fixed);
  o.precision(4);
  o << "%!PS-Adobe-3.0 EPSF-3.0\n";
  o << "%%Creator: schottky figure_export\n";
  if (!F.title.empty()) o << "%%Title: " << F.title << "\n";
  o << "%%BoundingBox: 0 0 " << (int)std::ceil(Wpt) << " " << (int)std::ceil(Hpt) << "\n";
  o << "%%HiResBoundingBox: 0 0 " << Wpt << " " << Hpt << "\n";
  o << "%%EndComments\n";
  o << "%%BeginProlog\n/M { moveto } bind def\n/L { lineto } bind def\n";
  o << "1 setlinejoin 1 setlinecap\n%%EndProlog\n";

  //background
  o << opt.bg_r/255.0 << " " << opt.bg_g/255.0 << " " << opt.bg_b/255.0
    << " setrgbcolor\n0 0 " << Wpt << " " << Hpt << " rectfill\n";

  if (!img.empty() && opt.embed_raster) {
    o << "gsave\n";
    o << "/picstr " << 3*img.w << " string def\n";
    o << "0 0 translate " << Wpt << " " << Hpt << " scale\n";
    o << img.w << " " << img.h << " 8 [" << img.w << " 0 0 " << -img.h
      << " 0 " << img.h << "]\n{currentfile picstr readhexstring pop} false 3 colorimage\n";
    //rows top-down, as the matrix above expects
    static const char* hex = "0123456789abcdef";
    std::string line;
    for (int y = img.h - 1; y >= 0; --y) {
      for (int x = 0; x < img.w; ++x) {
        unsigned char r, g, b;
        img.get_pixel(x, y, &r, &g, &b);
        unsigned char v[3] = {r, g, b};
        for (int c = 0; c < 3; ++c) {
          line += hex[(v[c] >> 4) & 0xf];
          line += hex[v[c] & 0xf];
        }
        if (line.size() >= 78) { o << line << "\n"; line.clear(); }
      }
    }
    if (!line.empty()) o << line << "\n";
    o << "grestore\n";
  }

  ps_emit_overlays(o, F, opt, T);
  o << "showpage\n%%EOF\n";

  std::ofstream fp(path.c_str(), std::ios::binary);
  if (!fp) { if (err) *err = "cannot open " + path + " for writing"; return false; }
  std::string s = o.str();
  fp.write(s.data(), (std::streamsize)s.size());
  if (!fp) { if (err) *err = "write failed on " + path; return false; }
  return true;
}

/*================================================================== 7. PDF ==*/

std::string pdf_escape(const std::string& s) { return ps_escape(s); }  //same rules

//a circle as four Beziers
void pdf_circle(std::ostringstream& o, double cx, double cy, double r) {
  const double K = 0.5522847498307933;
  o << cx + r << " " << cy << " m\n";
  o << cx + r << " " << cy + K*r << " " << cx + K*r << " " << cy + r << " " << cx << " " << cy + r << " c\n";
  o << cx - K*r << " " << cy + r << " " << cx - r << " " << cy + K*r << " " << cx - r << " " << cy << " c\n";
  o << cx - r << " " << cy - K*r << " " << cx - K*r << " " << cy - r << " " << cx << " " << cy - r << " c\n";
  o << cx + K*r << " " << cy - r << " " << cx + r << " " << cy - K*r << " " << cx + r << " " << cy << " c\n";
  o << "h\n";
}

void pdf_style(std::ostringstream& o, const Style& st, double lw_default) {
  double lw = (st.line_width > 0 ? st.line_width : lw_default);
  o << st.r << " " << st.g << " " << st.b << " RG\n";
  o << st.r << " " << st.g << " " << st.b << " rg\n";
  o << lw << " w\n";
  if (st.dashed) o << "[" << 3*lw << " " << 2*lw << "] 0 d\n";
  else           o << "[] 0 d\n";
}

bool write_pdf(const Figure& F, const Options& opt, const Raster& img,
               double Wpt, double Hpt, const std::string& path, std::string* err) {
  Xform T;
  T.x0 = F.x0; T.y0 = F.y0;
  T.sx = Wpt/(F.x1 - F.x0);
  T.sy = Hpt/(F.y1 - F.y0);

  bool have_img = (!img.empty() && opt.embed_raster);
  bool have_text = !F.texts.empty();

  /*--- the content stream ---*/
  std::ostringstream c;
  c.setf(std::ios::fixed);
  c.precision(4);
  c << opt.bg_r/255.0 << " " << opt.bg_g/255.0 << " " << opt.bg_b/255.0 << " rg\n";
  c << "0 0 " << Wpt << " " << Hpt << " re f\n";
  if (have_img)
    c << "q " << Wpt << " 0 0 " << Hpt << " 0 0 cm /Im0 Do Q\n";

  for (size_t k = 0; k < F.paths.size(); ++k) {
    const Path& p = F.paths[k];
    if (p.pts.size() < 2) continue;
    pdf_style(c, p.st, opt.line_width_pt);
    c << T.to_x(p.pts[0].x) << " " << T.to_y(p.pts[0].y) << " m\n";
    for (size_t m = 1; m < p.pts.size(); ++m)
      c << T.to_x(p.pts[m].x) << " " << T.to_y(p.pts[m].y) << " l\n";
    if (p.closed) c << "h\n";
    c << (p.st.filled ? "f\n" : "S\n");
  }
  for (size_t k = 0; k < F.disks.size(); ++k) {
    const Disk& d = F.disks[k];
    pdf_style(c, d.st, opt.line_width_pt);
    pdf_circle(c, T.to_x(d.c.x), T.to_y(d.c.y), d.r*T.sx);
    c << (d.st.filled ? "f\n" : "S\n");
  }
  for (size_t k = 0; k < F.dots.size(); ++k) {
    const Dot& d = F.dots[k];
    pdf_style(c, d.st, opt.line_width_pt);
    pdf_circle(c, T.to_x(d.c.x), T.to_y(d.c.y), d.r_pt);
    c << "f\n";
  }
  for (size_t k = 0; k < F.texts.size(); ++k) {
    const Text& t = F.texts[k];
    //PDF has no stringwidth operator, so anchor with the real Helvetica metrics
    double adv = helvetica_width(t.s, t.size_pt);
    double x = T.to_x(t.p.x);
    if (t.anchor == 1) x -= adv/2.0;
    else if (t.anchor == 2) x -= adv;
    c << t.st.r << " " << t.st.g << " " << t.st.b << " rg\n";
    c << "BT /F1 " << t.size_pt << " Tf " << x << " " << T.to_y(t.p.y) << " Td ("
      << pdf_escape(t.s) << ") Tj ET\n";
  }
  std::string content = c.str();

  /*--- objects.  1 catalog, 2 pages, 3 page, 4 contents, then image, then font ---*/
  int obj_img  = have_img  ? 5 : 0;
  int obj_font = have_text ? (have_img ? 6 : 5) : 0;
  int nobj = 4 + (have_img ? 1 : 0) + (have_text ? 1 : 0);

  std::vector<unsigned char> imgdata;
  if (have_img) {
    std::vector<unsigned char> raw;
    raw.reserve((size_t)img.w*img.h*3);
    for (int y = img.h - 1; y >= 0; --y)               //PDF images run top-down
      for (int x = 0; x < img.w; ++x) {
        unsigned char r, g, b;
        img.get_pixel(x, y, &r, &g, &b);
        raw.push_back(r); raw.push_back(g); raw.push_back(b);
      }
    zlib_wrap(raw.empty() ? (const unsigned char*)"" : &raw[0], raw.size(), imgdata);
  }

  std::string out;
  std::vector<size_t> off(nobj + 1, 0);
  out += "%PDF-1.4\n";
  out += "%\xe2\xe3\xcf\xd3\n";                        //binary marker

  char buf[512];

  off[1] = out.size();
  out += "1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n";

  off[2] = out.size();
  out += "2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n";

  off[3] = out.size();
  {
    std::ostringstream r;
    r.setf(std::ios::fixed); r.precision(4);
    r << "3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 " << Wpt << " " << Hpt
      << "] /Contents 4 0 R /Resources << ";
    if (have_img)  r << "/XObject << /Im0 " << obj_img << " 0 R >> ";
    if (have_text) r << "/Font << /F1 " << obj_font << " 0 R >> ";
    r << ">> >>\nendobj\n";
    out += r.str();
  }

  off[4] = out.size();
  /* snprintf, not sprintf: the values here are all short numbers so nothing could
     actually overflow buf, but sprintf is deprecated on macOS and warned on every build,
     which is noise a reader has to learn to ignore. */
  std::snprintf(buf, sizeof buf, "4 0 obj\n<< /Length %lu >>\nstream\n", (unsigned long)content.size());
  out += buf;
  out += content;
  out += "endstream\nendobj\n";

  if (have_img) {
    off[obj_img] = out.size();
    std::snprintf(buf, sizeof buf,
      "%d 0 obj\n<< /Type /XObject /Subtype /Image /Width %d /Height %d "
      "/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length %lu >>\nstream\n",
      obj_img, img.w, img.h, (unsigned long)imgdata.size());
    out += buf;
    out.append((const char*)(imgdata.empty() ? (const unsigned char*)"" : &imgdata[0]),
               imgdata.size());
    out += "\nendstream\nendobj\n";
  }

  if (have_text) {
    off[obj_font] = out.size();
    std::snprintf(buf, sizeof buf, "%d 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica "
                      "/Encoding /WinAnsiEncoding >>\nendobj\n", obj_font);
    out += buf;
  }

  size_t xref = out.size();
  std::snprintf(buf, sizeof buf, "xref\n0 %d\n", nobj + 1);
  out += buf;
  out += "0000000000 65535 f \n";
  for (int i = 1; i <= nobj; ++i) {
    std::snprintf(buf, sizeof buf, "%010lu 00000 n \n", (unsigned long)off[i]);
    out += buf;
  }
  std::snprintf(buf, sizeof buf, "trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%lu\n%%%%EOF\n",
               nobj + 1, (unsigned long)xref);
  out += buf;

  std::ofstream fp(path.c_str(), std::ios::binary);
  if (!fp) { if (err) *err = "cannot open " + path + " for writing"; return false; }
  fp.write(out.data(), (std::streamsize)out.size());
  if (!fp) { if (err) *err = "write failed on " + path; return false; }
  return true;
}

} //anonymous namespace

/*=================================================== 8. the public surface ==*/

Style::Style() {
  r = g = b = 0.0;
  line_width = 0.0;
  filled = false;
  dashed = false;
}

Style Style::stroke(double R, double G, double B, double lw) {
  Style s; s.r = R; s.g = G; s.b = B; s.line_width = lw; return s;
}

Style Style::fill(double R, double G, double B) {
  Style s; s.r = R; s.g = G; s.b = B; s.filled = true; return s;
}

Style Style::dash(double R, double G, double B, double lw) {
  Style s; s.r = R; s.g = G; s.b = B; s.line_width = lw; s.dashed = true; return s;
}

void Raster::set_size(int W, int H, unsigned char fill) {
  w = W; h = H;
  rgb.assign((size_t)3*W*H, fill);
}

Figure::Figure() { x0 = y0 = 0.0; x1 = y1 = 1.0; }

void Figure::set_window(double X0, double Y0, double X1, double Y1) {
  x0 = X0; y0 = Y0; x1 = X1; y1 = Y1;
  //a degenerate window would divide by zero downstream
  if (x1 - x0 <= 0.0) x1 = x0 + 1e-12;
  if (y1 - y0 <= 0.0) y1 = y0 + 1e-12;
}

void Figure::add_path(const std::vector<Point2d<double> >& p, const Style& s, bool closed) {
  Path q; q.pts = p; q.st = s; q.closed = closed;
  paths.push_back(q);
}

void Figure::add_segment(double ax, double ay, double bx, double by, const Style& s) {
  Path q;
  q.pts.push_back(Point2d<double>(ax, ay));
  q.pts.push_back(Point2d<double>(bx, by));
  q.st = s;
  paths.push_back(q);
}

void Figure::add_disk(double cx, double cy, double r, const Style& s) {
  Disk d; d.c = Point2d<double>(cx, cy); d.r = r; d.st = s;
  disks.push_back(d);
}

void Figure::add_dot(double cx, double cy, const Style& s, double r_pt) {
  Dot d; d.c = Point2d<double>(cx, cy); d.r_pt = r_pt; d.st = s;
  dots.push_back(d);
}

void Figure::add_text(double px, double py, const std::string& s, const Style& st,
                      double size_pt, int anchor) {
  Text t; t.p = Point2d<double>(px, py); t.s = s; t.st = st;
  t.size_pt = size_pt; t.anchor = anchor;
  texts.push_back(t);
}

void Figure::add_circle(double cx, double cy, double r, const Style& s, int nseg) {
  if (nseg < 8) nseg = 8;
  std::vector<Point2d<double> > p;
  circle_pts(cx, cy, r, nseg, p);
  add_path(p, s, true);
}

Options::Options() {
  format = PNG;
  width_pt = 360.0;              //5 inches
  raster_px = 0;                 //keep the raster's own resolution
  vector_overlays = true;
  transparent = false;
  bg_r = bg_g = bg_b = 255;
  draw_frame = true;
  draw_axes = false;
  line_width_pt = 0.6;
  embed_raster = true;
}

const char* format_name(Options::Format f) {
  switch (f) {
    case Options::PNG: return "PNG";
    case Options::EPS: return "EPS";
    case Options::PDF: return "PDF";
  }
  return "?";
}

bool format_from_extension(const std::string& path, Options::Format* out) {
  size_t d = path.rfind('.');
  if (d == std::string::npos) return false;
  std::string e = path.substr(d);
  for (size_t i = 0; i < e.size(); ++i)
    if (e[i] >= 'A' && e[i] <= 'Z') e[i] = (char)(e[i] - 'A' + 'a');
  if (e == ".png") { *out = Options::PNG; return true; }
  if (e == ".eps" || e == ".epsf" || e == ".ps") { *out = Options::EPS; return true; }
  if (e == ".pdf") { *out = Options::PDF; return true; }
  return false;
}

bool write_figure(const Figure& F_in, const Options& opt, const std::string& path,
                  std::string* err) {
  if (F_in.x1 <= F_in.x0 || F_in.y1 <= F_in.y0) {
    if (err) *err = "write_figure: degenerate window; call Figure::set_window";
    return false;
  }
  if (F_in.raster.empty() && !F_in.has_overlays()) {
    if (err) *err = "write_figure: nothing to draw (no raster, no overlays)";
    return false;
  }

  Figure F = F_in;
  add_decorations(F, opt);

  //aspect ratio is set by the mathematical window, never by the output size
  double aspect = (F.y1 - F.y0)/(F.x1 - F.x0);
  double Wpt = opt.width_pt > 0 ? opt.width_pt : 360.0;
  double Hpt = Wpt*aspect;

  //the raster we will actually emit
  Raster img;
  if (!F.raster.empty()) {
    if (opt.raster_px > 0 && opt.raster_px != F.raster.w) {
      int W = opt.raster_px;
      int H = (int)(W*aspect + 0.5);
      if (H < 1) H = 1;
      resample(F.raster, W, H, img);
    } else {
      img = F.raster;
    }
  }

  if (opt.format == Options::PNG) {
    //PNG is pure raster: make a canvas, then rasterise the overlays onto it
    int W, H;
    if (!img.empty()) { W = img.w; H = img.h; }
    else {
      W = opt.raster_px > 0 ? opt.raster_px : 1000;
      H = (int)(W*aspect + 0.5);
      if (H < 1) H = 1;
    }
    Raster canvas;
    if (!img.empty() && img.w == W && img.h == H) canvas = img;
    else if (!img.empty()) resample(img, W, H, canvas);
    else canvas.set_size(W, H, 255);
    if (img.empty()) {
      for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i) canvas.set_pixel(i, j, opt.bg_r, opt.bg_g, opt.bg_b);
    }
    Canvas C;
    C.R = &canvas;
    C.T.x0 = F.x0; C.T.y0 = F.y0;
    C.T.sx = W/(F.x1 - F.x0);
    C.T.sy = H/(F.y1 - F.y0);
    C.px_per_pt = W/Wpt;
    rasterise_overlays(F, opt, C);
    //transparency for a truecolor PNG means "this one color is see-through",
    //so it can only be honoured where there is a background color to name
    return write_png_impl(path, W, H, canvas.rgb, opt.transparent,
                          opt.bg_r, opt.bg_g, opt.bg_b, err);
  }

  //EPS and PDF.  If the caller asked for rasterised overlays, burn them into the
  //image first and then emit a Figure with no overlays.
  if (!opt.vector_overlays && !img.empty()) {
    Canvas C;
    C.R = &img;
    C.T.x0 = F.x0; C.T.y0 = F.y0;
    C.T.sx = img.w/(F.x1 - F.x0);
    C.T.sy = img.h/(F.y1 - F.y0);
    C.px_per_pt = img.w/Wpt;
    rasterise_overlays(F, opt, C);
    Figure G;
    G.set_window(F.x0, F.y0, F.x1, F.y1);
    G.title = F.title;
    if (opt.format == Options::EPS) return write_eps(G, opt, img, Wpt, Hpt, path, err);
    return write_pdf(G, opt, img, Wpt, Hpt, path, err);
  }

  if (opt.format == Options::EPS) return write_eps(F, opt, img, Wpt, Hpt, path, err);
  return write_pdf(F, opt, img, Wpt, Hpt, path, err);
}

bool write_auto(const Figure& F, Options opt, const std::string& path, std::string* err) {
  Options::Format f;
  if (!format_from_extension(path, &f)) {
    if (err) *err = "cannot tell the format from the name '" + path +
                    "'; use .png, .eps or .pdf";
    return false;
  }
  opt.format = f;
  return write_figure(F, opt, path, err);
}

void viridis(double t, unsigned char rgb[3]) {
  /* Keep in step with RAMP in render_funddom.py. */
  static const double S[11][3] = {
    {0.267004, 0.004874, 0.329415},
    {0.282623, 0.140926, 0.457517},
    {0.253935, 0.265254, 0.529983},
    {0.206756, 0.371758, 0.553117},
    {0.163625, 0.471133, 0.558148},
    {0.127568, 0.566949, 0.550556},
    {0.134692, 0.658636, 0.517649},
    {0.266941, 0.748751, 0.440573},
    {0.477504, 0.821444, 0.318195},
    {0.741388, 0.873449, 0.149561},
    {0.993248, 0.906157, 0.143936}
  };
  if (!(t > 0.0)) t = 0.0;          /* also catches nan */
  if (t > 1.0) t = 1.0;
  double x = t*10.0;
  int i = (int)x;
  if (i > 9) i = 9;
  double u = x - i;
  for (int k=0; k<3; ++k) {
    double v = S[i][k] + u*(S[i+1][k] - S[i][k]);
    int q = (int)(255.0*v + 0.5);
    if (q < 0) q = 0;
    if (q > 255) q = 255;
    rgb[k] = (unsigned char)q;
  }
}

} //namespace figexp
