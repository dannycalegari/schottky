/* connected.cc	

	Implements circle algorithm to test for connectedness.

*/

//the input disks are c1, a1, = (center of the disk, vector pointing to 
bool ifs::old_circles_intersect(cpx c1, cpx a1, cpx c2, cpx a2, double R, int d){
	cpx c1p, a1p, c1q, a1q, c2p, a2p, c2q, a2q;
	
	//	draw_circle(cpx_to_point(c1),cpx_to_radius(abs(a1)*R),0xFF0000);
	//	draw_circle(cpx_to_point(c2),cpx_to_radius(abs(a2)*R),0xFF0000);


//	cout << "testing at depth " << d << "\n";
	if(abs(c2-c1) > abs(a1)*R+abs(a2)*R){	// are circles naively disjoint?
		return(false);
	} else {
		if(d>0){
			// need to subdivide
			
			/* this shortcut allows fast calculation
			of connectedness locus but not internal structure
			of disconnectedness locus (and not much fine structure) */
			
		//	if(exit_depth>depth){
		//		return(true);
		//	};
						
			exit_depth++;
			// 4 new circles
			c1p = c1-a1/2.0;
			a1p = a1*z;
			c1q = c1+a1/2.0;
			a1q = a1*w;
			c2p = c2-a2/2.0;
			a2p = a2*z;
			c2q = c2+a2/2.0;
			a2q = a2*w;
			
			// have 4 chances for an intersection 
			// c1p,a1p,c2p,a2p ; c1p,a1p,c2q,a2q ; c1q,a1q,c2p,a2p ; c1q,a1q,c2q,a2q
			if(old_circles_intersect(c1p,a1p,c2p,a2p,R,d-1) || old_circles_intersect(c1p,a1p,c2q,a2q,R,d-1) ||
				old_circles_intersect(c1q,a1q,c2p,a2p,R,d-1) || old_circles_intersect(c1q,a1q,c2q,a2q,R,d-1)){
				// if 1 of the 4 pairs recursively intersect after refinement to maximum depth
				return(true);
			} else {
				return(false);
			};
		} else {
			return(true);
		};
	};
};

//recursive test if the circles intersect
//the disks are centered at the centers, and you get to the image of the
//center by adding the z_img_i w_img_i vectors
bool ifs::circles_intersect(cpx center_1, cpx z_img_1, cpx w_img_1, double R1,
                            cpx center_2, cpx z_img_2, cpx w_img_2, double R2, int d) {
  if ( abs(center_1 - center_2) > R1 + R2) {
    return false;
  } else if (d<=0) {
    return true;
  }
  cpx c1z, c1w, zi1z, zi1w, wi1z, wi1w;
  cpx c2z, c2w, zi2z, zi2w, wi2z, wi2w;
  double R1z, R1w, R2z, R2w;
  c1z = center_1 + z_img_1;
  zi1z = z_img_1*z;
  wi1z = w_img_1*z;
  R1z = R1*az;
  
  c1w = center_1 + w_img_1;
  zi1w = z_img_1*w;
  wi1w = w_img_1*w;
  R1w = R1*aw;
  
  c2z = center_2 + z_img_2;
  zi2z = z_img_2*z;
  wi2z = w_img_2*z;
  R2z = R2*az;
  
  c2w = center_2 + w_img_2;
  zi2w = z_img_2*w;
  wi2w = w_img_2*w;
  R2w = R2*aw;
  
  ++exit_depth;
  
  //there's no point in checking 
  //the images inside a single disk, because those will always be connected
  //to a higher level than we're checking here, so we check all 4 possible
  //pairs across the disks
  if (circles_intersect( c1z, zi1z, wi1z, R1z, c2z, zi2z, wi2z, R2z, d-1 ) ||
      circles_intersect( c1z, zi1z, wi1z, R1z, c2w, zi2w, wi2w, R2w, d-1 ) ||
      circles_intersect( c1w, zi1w, wi1w, R1w, c2z, zi2z, wi2z, R2z, d-1 ) ||
      circles_intersect( c1w, zi1w, wi1w, R1w, c2w, zi2w, wi2w, R2w, d-1 ) ){
    return true;
  }
  return false;
}


bool ifs::circ_connected(double r){ 	// circle algorithm to test for connectedness
	double R;
	if (r<0) {
	  R=max( abs(z-1.0)/(2.0*(1.0-az)), abs(w-1.0)/(2.0*(1.0-aw)) );
	} else {
	  R = r;
	}
	
//	draw_circle(cpx_to_point(0.5),cpx_to_radius(R),0xFF0000);
        //std::cout << "Checking connectedness with radius " << R << "\n";
	exit_depth=0;	// initialize exit_depth
	bool ans = circles_intersect(0.5*z, (z*z-z)/2.0, (z-z*w)/2.0, az*R,
                                   1.0-w/2.0, (z*w-w)/2.0, (w-w*w)/2.0, aw*R, depth);
        //std::cout << "Got " << ans << "\n";
	return ans;
}



bool ifs::is_connected(int d, int& difficulty) {
  if ((abs(z-w) < 0.000001) && az > 1.0/sqrt(2.0)) {
    difficulty = 0;
    return true;
  }
  double min_r;
  if (!minimal_enclosing_radius(min_r)) {
    difficulty = 1;
    return true;
  }
  //std::cout << "Checking if " << z << " " << w << " is connected with minimal radius " << min_r << "\n";
  int old_depth = depth;
  depth = d;
  bool ans = circ_connected(min_r);
  difficulty = exit_depth;
  depth = old_depth;
  return ans;
}



bool ifs::contains_point_recurse(const cpx& pt, const Ball& b, int d) {
  if (abs(pt - b.center) > b.radius) {
    return false;
  }
  if (d <= 0) return true;
  Ball bz = act_on_right(0, b);
  Ball bw = act_on_right(1, b);
  ++exit_depth;
  return (contains_point_recurse(pt, bz, d-1) || contains_point_recurse(pt, bw, d-1));
}
  





bool ifs::contains_point(cpx pt, double r) {
  double R;
  if (r<0) {
    R=max( abs(z-1.0)/(2.0*(1.0-az)), abs(w-1.0)/(2.0*(1.0-aw)) );
  } else {
    R = r;
  }
  exit_depth = 0;
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,R);
  bool ans = contains_point_recurse(pt, initial_ball, depth);
  return ans;
}



bool ifs::contains_half(int d, int& difficulty) {
  double min_r;
  if (!minimal_enclosing_radius(min_r)) {
    difficulty = 0;
    return true;
  }
  int old_depth = depth;
  depth = d;
  bool ans = contains_point(0.5, min_r);
  depth = old_depth;
  difficulty = exit_depth;
  return ans;
}







