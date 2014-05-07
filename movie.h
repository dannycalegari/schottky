#ifndef __MOVIE__
#define __MOVIE__

#include <vector>
#include <string>

#include "cpx.h"
#include "ifs.h"

#define PATH_TO_FFMPEG "../ffmpeg/"

bool write_bitmap(const std::vector<std::vector<Point3d<unsigned char> > >& p, 
                  const std::string& filename);
                  
bool ifs_movie_from_path(const ifs& IFS, 
                         const std::vector<cpx>& path, 
                         bool close_path,
                         const std::string& filename,
                         const cpx& region_ll,
                         const cpx& region_ur,
                         int depth,
                         bool with_mandlebrot,
                         cpx* mand_ll,
                         cpx* mand_ur,
                         std::vector<std::vector<bool> >* mand_connected_grid,
                         int pixel_w,
                         int pixel_h,
                         int fps, 
                         int seconds,
                         int verbose);



#endif