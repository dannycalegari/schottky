#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "point.h"

bool write_bitmap(const std::vector<std::vector<Point3d<char> > >& p, 
                  const std::string& filename) {
  
  if (p.size() == 0 || p[0].size() == 0) return false;
  
  FILE* f = NULL;
  std::string s = filename;
  f = fopen(s.c_str(), "wb");
  if (f == NULL) return false;
  
  int w = p.size(); //width in pixels
  int h = p[0].size(); //height in pixels
  int row_remainder = (3*w) % 4;
  int padding_bytes_per_row = (row_remainder == 0 ? 0 : 4-row_remainder);
  int row_btyes = 3*w + padding_bytes_per_row;
  int data_length = h*row_btyes;
  int total_file_length = data_length + 54;
  
  unsigned char* raw_data = (unsigned char*)malloc(data_length*sizeof(unsigned char));
  
  unsigned char BMP_header[14] = {'B','M',                 //bitmap code
                                  'X','X','X','X',          //will fill in total file size (bytes) here
                                   0, 0, 0, 0,             //unused
                                  54, 0, 0, 0 };           //offset where pixels start (14BMP header + 40DIB header)
  unsigned char DIB_header[40] = { 40, 0, 0, 0,            //bytes in DIB header
                                 'W','W','W','W',         //will fill in width in pixels
                                 'H','H','H','H',         //will fill in height in pixels (positive for bottom-to-top order)
                                  1 , 0,                  //one color plane used?
                                  24, 0,                  //bits per pixels (BGR order!)
                                  0, 0, 0, 0,             //no compressiong
                                 'D','D','D','D',         //will fill in with size (bytes) of bitmap data (including padding)
                                 0x13, 0x0B,0,0,          //2835 pixels/meter horizontal
                                 0x13, 0x0B,0,0,          //2835 pixels/meter vertical
                                 0,0,0,0,                 //no color palette
                                 0,0,0,0 };               //no important colors
  
  //write in the correct total file size
  BMP_header[2] = (unsigned char)total_file_length;
  BMP_header[3] = (unsigned char)(total_file_length >> 8);
  BMP_header[4] = (unsigned char)(total_file_length >> 16);
  BMP_header[5] = (unsigned char)(total_file_length >> 24);
  
  //write in the width and height in pixels
  DIB_header[4] = (unsigned char)w;
  DIB_header[5] = (unsigned char)(w >> 8);
  DIB_header[6] = (unsigned char)(w >> 16);
  DIB_header[7] = (unsigned char)(w >> 24);
  DIB_header[8] = (unsigned char)h;
  DIB_header[9] = (unsigned char)(h >> 8);
  DIB_header[10] = (unsigned char)(h >> 16);
  DIB_header[11] = (unsigned char)(h >> 24);
  
  //write in the data size (including padding)
  DIB_header[20] = (unsigned char)data_length;
  DIB_header[21] = (unsigned char)(data_length >> 8);
  DIB_header[22] = (unsigned char)(data_length >> 16);
  DIB_header[23] = (unsigned char)(data_length >> 24);
  
  //create the bmp data array
  int current_data_index = 0;
  for (int i=0; i<h; ++i) {
    for (int j=0; j<w; ++j) {
      raw_data[current_data_index] = p[j][i].z;
      raw_data[current_data_index+1] = p[j][i].y;
      raw_data[current_data_index+2] = p[j][i].x;
      current_data_index += 3;
    }
    for (int j=0; j<padding_bytes_per_row; ++j) {
      raw_data[current_data_index] = 0;
      current_data_index++;
    }
  }
  
  //write the file
  fwrite(BMP_header, sizeof(char), 14, f);
  fwrite(DIB_header, sizeof(char), 40, f);
  fwrite(raw_data, sizeof(char), data_length, f);
  fclose(f);
  
  return 0;
} 







int main(int argc, char* argv[]) {
  std::stringstream T;
  T.str("");
  std::vector<std::vector<Point3d<char> > > arr(100);
  for (int i=0; i<100; ++i) {
    arr[i].resize(100);
    for (int j=0; j<100; ++j) {
      arr[i][j] = Point3d<char>(0,0,0);
    }
  }
  for (int i=0; i<100; ++i) {
    if (i > 0) arr[i-1][i-1] = Point3d<char>(0,0,0);
    arr[i][i] = Point3d<char>(255,255,255);
    T.str("");
    T << "bitmap_test" << i << ".bmp";
    std::string s = T.str();
    write_bitmap(arr, s);
  }
  return 0;
}



























