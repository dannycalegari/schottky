/* ifs.cc */

class ifs{
	public:
		// data defining IFS and display and adjust options
		/* IFS is generated by two similarities; one scales by z and is centered at 0,
		the other scales by w and is centered at 1. */
		cpx z,w;		// parameters for IFS
		bool color_ifs;	// draw fL and gL in different colors
		int sync;		// sync=0 w is arbitrary, sync=1 w=z, sync=2 w=-z
		int depth;		// depth to recurse to draw or detect connectedness
		int exit_depth;	// size of tree to detect connectedness
		cpx seed;		// initial seed for IFS in ifs mode
		double step;	// step size to adjust z or w in ifs mode
		cpx center;		// center of screen in mandelbrot mode
		double wind;	// size of window in mandelbrot mode
		int mesh;		// size of mesh in mandelbrot mode in pixels
		int mode;		// draw mode: 0 for limit set, 1 for mandelbrot set
		
		bool disconnection_depth; // whether to draw the depth for disconnected sets
						
		// default initial values
		
		void initialize(cpx a, cpx b);
		
		// apply generator; should really make this function more flexible
		
		cpx iterate(int index, cpx u);
				
		// test connectedness of limit set
		
		bool circles_intersect(cpx c1, cpx a1, cpx c2, cpx a2, double R, int d); // recursive test; do circles intersect?
		bool circ_connected();	// circle algorithm to test for connectedness
				
		// draw limit set
		
		void draw_dots(int depth, cpx u);
		void draw_color_dots(int d, cpx u, long color);

		void draw_limit_set();
		
		// draw mandelbrot set
		
		void zoom(point p);
		void draw_mandelbrot_set();
		
		// overall draw
		
		void draw();
		
		// user interface
			
		void user_interface();
};


void ifs::initialize(cpx a, cpx b){
	// initialize z to a and w to b
	z=a;
	w=b;
	sync=0;
	color_ifs=true;
	disconnection_depth=false;
	step=0.01;	// size of adjustments to z and w
	seed=0.0;	// initial seed point for IFS
	center=0.0;	// in mandelbrot mode; center of screen, size of window, and mesh of accuracy
	wind=1.0;
	mesh=1;
	depth=12;	// depth to iterate IFS or detect connectedness to
};

cpx ifs::iterate(int index, cpx u){
	// apply generator to u
	if(index==0){
		return(u*z);
	} else {
		return(((u-1.0)*w)+1.0);
	};
};

void ifs::draw_dots(int d, cpx u){
	// draw IFS to depth d with seed u
	if(d<=0){
		draw_point(cpx_to_point(u),0);
	} else {
		draw_dots(d-1, iterate(0,u));
		draw_dots(d-1, iterate(1,u));
	};
};

void ifs::draw_color_dots(int d, cpx u, long color){
	// draw IFS to depth d with seed u
	if(d<=0){
		draw_point(cpx_to_point(u), color);
	} else {
		draw_color_dots(d-1, iterate(0,u), 0xFFAA00);
		draw_color_dots(d-1, iterate(1,u), 0x00AAFF);
	};
};

void ifs::draw_limit_set(){
	// draw IFS and write values of parameters to screen
	stringstream T;
	point p;
	
	if(sync==1){
		w=z;
	} else if(sync==2){
		w=conj(z);
	};
	erase_field();
	p.x=1200;
	p.y=20;
	T << "limit set (IFS attractor)";
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;	
	T << "z (adjust with [arrows]): " << z;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;	
	T << "w (adjust with [ijkl]): " << w;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "depth (adjust with [d/e]): " << depth;
	draw_text(p,T,0x000000);
	T.str("");	
	p.y=p.y+20;
	T << "step (adjust with [a/s]): " << step;
	draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y+20;
	T << "sync w free, w=z, w=conj(z) (cycle with [v]): " << sync;
	draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y+20;
	T << "2-color IFS (toggle with [c]): ";
	if(color_ifs){
		T << "on";
	} else {
		T << "off";
	};
	draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y+20;
	T << "connected: ";
	if(circ_connected()){
		T << "yes";
	} else {
		T << "no";
	};
	draw_text(p,T,0x000000);
	T.str("");	
	p.y=p.y+20;
	T << "toggle ifs/mandelbrot with [b]";
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "quit with [q]";
	draw_text(p,T,0x000000);
	T.str("");	
	if(color_ifs){
		draw_color_dots(depth,seed,0x000000);	// actual IFS is drawn
	} else {
		draw_dots(depth,seed);
	};
	draw_dot(cpx_to_point(z),0xFF0000);	// z
	draw_dot(cpx_to_point(1.0-w),0x00FF00);	// 1-w
	draw_dot(cpx_to_point(0.0),0xBDB76B);	// 0
	draw_dot(cpx_to_point(1.0),0xBDB76B);	// 1
};

void ifs::zoom(point p){
	// in mandelbrot
	cpx cc;
	cc=point_to_cpx(p);
	center=(cc*wind)+center;
	wind=wind/2.0;
	mesh=mesh*2;
};

void ifs::draw_mandelbrot_set(){
	// draw mandelbrot set (connectedness locus) and write values of parameters to screen
	stringstream T;
	point p,q;
	cpx zz;
	int i,j;
	double x,y;
	
	erase_field();
	
	p.x=1200;	// 1200
	p.y=20;
	T << "connectedness locus";
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "slice: ";
	if(sync==0){
		T << "w=" << w << "\n";
	} else if(sync==1){
		T << "w=z\n";
	} else if(sync==2){
		T << "w=conj(z)\n";
	};
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "cycle slice type with [v]\n";
	draw_text(p,T,0x000000);
	T.str("");		
	p.y=p.y+20;
	T << "depth (adjust with [d/e]): " << depth;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "center (adjust by mouse click): " << center;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "window (zoom out with [g]): " << wind;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "mesh (adjust with [n/m]): " << mesh;
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "toggle ifs/mandelbrot with [b]";
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "toggle disconnection depth with [f]";
	draw_text(p,T,0x000000);
	T.str("");
	p.y=p.y+20;
	T << "quit with [q]";
	draw_text(p,T,0x000000);
	T.str("");	
	
	for(i=0;i<1024;i=i+mesh){
		for(j=0;j<1024;j=j+mesh){
			q.x=i;
			q.y=j;
			x=2.0*wind*((double) i)/1024.0;
			y=2.0*wind*((double) 1024-j)/1024.0;
			
			z=center-wind-(wind*I)+x+(y*I);
			if(sync==1){
				w=z;	// diagonal slice; could define other slices (eg w=constant,w=conj(z))
			} else if(sync==2){
				w=conj(z);
			};
			if(abs(z)>1.0){	// could truncate this to sqrt(1/2) actually
				draw_box(q,mesh,0x000000);
			} else { // if(abs(z)>0.5){
				if(circ_connected()){
					draw_box(q,mesh,0x000001*exit_depth);
				} else if (disconnection_depth) {
					draw_box(q,mesh,0x010000*exit_depth);
				}
			};
			
	
		};
	};
};

void ifs::draw(){
	if(mode==0){
		draw_limit_set();
	} else if(mode==1){
		draw_mandelbrot_set();
	};
};