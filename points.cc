/* points.cc definitions and basic functions of points	*/

point cpx_to_point(cpx w){
    point p;
	if(SMALLSCREEN){
	    p.x = 256 + (int) 256.0*w.real();
	    p.y = 256 - (int) 256.0*w.imag();
	} else {
		p.x = 512 + (int) 512.0*w.real(); // 512
		p.y = 512 - (int) 512.0*w.imag(); // 512
	};
    return(p);
};

int cpx_to_radius(cpx w){
	int r;
	r= (int) 512.0*abs(w);
    return(r);
};


cpx point_to_cpx(point p){
	cpx w;
	w.real() = ((double) (p.x-512))/512.0;
	w.imag() = ((double) (p.y-512))/-512.0;
	return(w);
};