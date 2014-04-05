
//draw all the dots at depth d with seed u (recursively)
void ifs::draw_dots(int d, cpx u){
	if(d<=0){
		X.draw_point(cpx_to_point(u),0);
	} else {
		draw_dots(d-1, iterate(0,u));
		draw_dots(d-1, iterate(1,u));
	}
}

//same as above, but colors for z and w
void ifs::draw_color_dots(int d, cpx u, long color){
	// draw IFS to depth d with seed u
	if(d<=0){
		X.draw_point(cpx_to_point(u), color);
	} else {
		draw_color_dots(d-1, iterate(0,u), 0xFFAA00);
		draw_color_dots(d-1, iterate(1,u), 0x00AAFF);
	}
}

//same as above, but draw them chunky starting at a given radius
void ifs::draw_color_chunky_dots(int d, cpx u, long color, double radius) {
	if (d<=0) {
	        //the radius should get multiplied by 512
	        double draw_radius = radius*drawing_radius;
	        if (draw_radius < 1) draw_radius = 1;
	        X.draw_disk(cpx_to_point(u), draw_radius, color);
	        //cout << "radius: " << radius;
	        //cout << "abs(z): " << abs(z);
	} else {
		draw_color_chunky_dots(d-1, iterate(0,u), 0xFFAA00, radius*az);
		draw_color_chunky_dots(d-1, iterate(1,u), 0x00AAFF, radius*aw);
	}
}

//draw the IFS and print out the data and stuff
void ifs::draw_limit_set(){
	std::stringstream T;
	Point2d<int> p;
	
	if(sync==1){
		w=z;
		aw = az;
	} else if(sync==2) {
		w=conj(z);
		aw = az;
	}
	X.erase_field();
	p.x=drawing_width+100;
	p.y=drawing_width - 40;
	T << "limit set (IFS attractor)";
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;	
	T << "z (adjust with [arrows]): " << z;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;	
	T << "w (adjust with [ijkl]): " << w;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "depth (adjust with [d/e]): " << depth;
	X.draw_text(p,T,0x000000);
	T.str("");	
	p.y=p.y-20;
	T << "step (adjust with [a/s]): " << step;
	X.draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y-20;
	T << "sync w free, w=z, w=conj(z) (cycle with [v]): " << sync;
	X.draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y-20;
	if (draw_trap_mode) {
		if (find_trap()) {
			T << "trap found (toggle trap mode with [t])";
			draw_trap();
		} else {
			T << "trap not found (toggle trap mode with [t])";
		}
	} else {
		T << "trap drawing disabled (toggle trap mode with [t])";
	}
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "2-color IFS (toggle with [c]): ";
	if(color_ifs){
		T << "on";
	} else {
		T << "off";
	};
	X.draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y-20;
	T << "chunky IFS (toggle with [x]): ";
	if (chunky_ifs) {
		T << "on (radius: " << minimal_enclosing_radius() << ")";
	} else {
		T << "off";
	}
	X.draw_text(p,T,0);
	T.str("");
	p.y=p.y-20;	
	T << "connected: ";
	if(circ_connected()){
		T << "yes";
	} else {
		T << "no";
	};
	X.draw_text(p,T,0x000000);
	T.str("");	
	p.y=p.y-20;
	T << "toggle ifs/mandelbrot with [b]";
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "quit with [q]";
	X.draw_text(p,T,0x000000);
	T.str("");	
	if(color_ifs){
		if (chunky_ifs) {
			double r = minimal_enclosing_radius();
			//cout << "Chunky radius: " << chunky_radius << "\n";
			draw_color_chunky_dots(depth, seed, 0x000000, r);
		} else {
			draw_color_dots(depth,seed,0x000000);	// actual IFS is drawn
		}
	} else {
		draw_dots(depth, seed);
	};
	X.draw_dot(cpx_to_point(z),0xFF0000);	// z
	X.draw_dot(cpx_to_point(1.0-w),0x00FF00);	// 1-w
	X.draw_dot(cpx_to_point(0.0),0xBDB76B);	// 0
	X.draw_dot(cpx_to_point(1.0),0xBDB76B);	// 1
};

void ifs::zoom(const Point2d<int>& p){
	// in mandelbrot
	cpx cc;
	cc=point_to_cpx(p);
	center=(cc*wind)+center;
	wind=wind/2.0;
	mesh=mesh*2;
};

void ifs::draw_mandelbrot_set(){
	// draw mandelbrot set (connectedness locus) and write values of parameters to screen
	std::stringstream T;
	Point2d<int> p,q;
	cpx I(0,1);
	cpx zz;
	int i,j;
	double x,y;
	
	X.erase_field();
	
	p.x=drawing_width + 100;
	p.y=drawing_width-40;
	T << "connectedness locus";
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "slice: ";
	if(sync==0){
		T << "w=" << w << "\n";
	} else if(sync==1){
		T << "w=z\n";
	} else if(sync==2){
		T << "w=conj(z)\n";
	};
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "cycle slice type with [v]\n";
	X.draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y-20;
	T << "depth (adjust with [d/e]): " << depth;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "center (adjust by mouse click): " << center;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "window (zoom out with [g]): " << wind;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "mesh (adjust with [n/m]): " << mesh;
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "toggle ifs/mandelbrot with [b]";
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "toggle disconnection depth with [f]";
	X.draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y-20;
	T << "quit with [q]";
	X.draw_text(p,T,0x000000);
	T.str("");	
	
	for(i=0; i<drawing_width; i=i+mesh){
		for(j=0; j<drawing_width; j=j+mesh){
			q.x=i;
			q.y=j;
			x=2.0*wind*double(i)/double(drawing_width);
			y=2.0*wind*double(j)/double(drawing_width);
			
			z=center-wind-(wind*I)+x+(y*I);
			if(sync==1){
				w=z;	// diagonal slice; could define other slices (eg w=constant,w=conj(z))
			} else if(sync==2){
				w=conj(z);
			};
			if(abs(z)>1.0){	// could truncate this to sqrt(1/2) actually
				X.draw_box(q,mesh,0x000000);
			} else { // if(abs(z)>0.5){
				if(circ_connected()){
					X.draw_box(q,mesh,0x000001*exit_depth);
				} else if (disconnection_depth) {
					X.draw_box(q,mesh,0x010000*exit_depth);
				}
			}
			
	
		}
	}
}

void ifs::draw(){
	if(mode==0){
		draw_limit_set();
	} else if(mode==1){
		draw_mandelbrot_set();
	}
}