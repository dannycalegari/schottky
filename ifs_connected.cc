/* connected.cc	

	Implements circle algorithm to test for connectedness.

*/

bool ifs::circles_intersect(cpx c1, cpx a1, cpx c2, cpx a2, double R, int d){
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
			if(circles_intersect(c1p,a1p,c2p,a2p,R,d-1) || circles_intersect(c1p,a1p,c2q,a2q,R,d-1) ||
				circles_intersect(c1q,a1q,c2p,a2p,R,d-1) || circles_intersect(c1q,a1q,c2q,a2q,R,d-1)){
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

bool ifs::circ_connected(double r){ 	// circle algorithm to test for connectedness
	double R;
	if (r<0) {
	  R=max(1.0/(2.0*(1.0-abs(z))),1.0/(2.0*(1.0-abs(w))));
	} else {
	  R = r;
	}
	
//	draw_circle(cpx_to_point(0.5),cpx_to_radius(R),0xFF0000);

	exit_depth=0;	// initialize exit_depth
	return(circles_intersect(0,z,1,w,R,depth));
};

