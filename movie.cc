#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>


#include "movie.h"



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
  
  return true;
} 


bool ifs_movie_from_path(const ifs& IFS, 
                         const std::vector<cpx>& path, 
                         bool close_path,
                         const std::string& filename,
                         const cpx& region_ll,
                         const cpx& region_ur,
                         int depth,
                         int pixel_w,
                         int pixel_h,
                         int fps, 
                         int seconds,
                         int verbose) {
  
  //create a copy of the ifs, because we're going to mess with it
  ifs IFS2(IFS.z, IFS.w, IFS.drawing_width, IFS.mode);
  
  //find the total path length and step size and stuff
  double total_path_len = 0;
  for (int i=0; i<(int)path.size()-1; ++i) {
    total_path_len += abs(path[i]-path[i+1]);
  }
  if (close_path) total_path_len += abs(path[path.size()-1] - path[0]);
  int num_frames = fps*seconds; 
  if (num_frames > 2000) {
    std::cout << "That's a lot of frames...\n";
    return false;
  }
  double length_per_frame = total_path_len / double(num_frames);
  if (verbose>0) {
    std::cout << "Box: " << region_ll << " " << region_ur << "\n";
    std::cout << "Path length: " << total_path_len << "\n";
    std::cout << "Total frames: " << num_frames << "\n";
    std::cout << "Distance per frame: " << length_per_frame << "\n";
  }
  
  //follow the path
  double length_done = 0;
  cpx z = path[0];
  int path_edge = 0;
  int frame_num = 0;
  std::stringstream T;
  std::vector<std::vector<Point3d<char> > > bmp(pixel_w, 
                                                std::vector<Point3d<char> >(pixel_h, 
                                                                            Point3d<char>(0,0,0)));
  do {
    //write out the bitmap
    if (verbose > 0) {
      std::cout << "\rWriting frame " << frame_num;
      std::cout.flush();
    }
    IFS2.set_params(z,z);
    IFS2.draw_ifs_to_array(bmp, region_ll, region_ur, depth);
    T.str("");
    T << filename << frame_num << ".bmp";
    write_bitmap(bmp, T.str());
    ++frame_num;
    
    
    //figure out the next z value
    double length_remaining = length_per_frame;
    while (true) {
      int next_path_ind = (path_edge == (int)path.size()-1 ? 0 : path_edge + 1);
      double to_next_endpoint = abs(z - path[next_path_ind]);
      if (to_next_endpoint > length_remaining) {
        cpx v = path[next_path_ind] - z;
        v = v/abs(v);
        z = z + length_remaining*v;
        break;
      }
      length_remaining -= to_next_endpoint;
      z = path[next_path_ind];
      path_edge = next_path_ind;
    }
    length_done += length_per_frame;
  } while (length_done < total_path_len);
  if (verbose>0) {
    std::cout << "\nWrote frame files\n";
    std::cout << "Calling movie encoder...";
    std::cout.flush();
  }
   
  //call ffmpeg to encode the movie
  std::stringstream command("");
  command << PATH_TO_FFMPEG << "ffmpeg -loglevel error -f image2 -r " << fps << " -i " << 
                                filename << "%d.bmp -c:v mpeg4 -qscale:v 8 " << 
                                filename << ".mp4";
  system(command.str().c_str());
  if (verbose>0) std::cout << "done\nErasing frame files...";
  
  //erase the frame files
  command.str("");
  command << "rm " << filename << "*.bmp";
  system(command.str().c_str());
  if (verbose>0) std::cout << "done\n";
  
  return true;
  
}

































