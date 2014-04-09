/* interface.cc */

void ifs::user_interface() {
	bool finished;
	Point2d<int> p,pp;
	cpx mz;
	std::stringstream T;
	int s;
	cpx I(0,1);
	
	finished=false;
	
	std::cout << "Entering interactive mode. Press [q] to exit.\n";
	std::cout << "(note: mouse must be focussed in X-window)\n";
	
	//initialize the graphics window
	//extra 580 for the right hand text
	X.initialize(drawing_width+350, drawing_width, 1, Point2d<float>(0,0));
	
	draw();
	
	while(finished==false){
	  XEvent report;
	  X.get_next_event(report);
		switch(report.type) {
			case ButtonPress:
				p=X.mouse_location();
				if(mode==1){	// mandelbrot mode
					zoom(p);
					draw();
				};	
				break;
			case MotionNotify:
				p=X.mouse_location();
				if(mode==1){	// mandelbrot mode
					mz=point_to_cpx(p);
					mz=(mz*wind)+center;
					p.x=drawing_width + 15;
					p.y=drawing_width/2;
					T.str("");
					T << "mouse location " << mz;
					pp=p;
					pp.y=pp.y-30;
					X.draw_filled_rectangle(pp, 300, 40, 0xFFFFFF);
					X.draw_text(p, T, 0x000000);
				};
				break;
			case KeyPress:
				if(XLookupKeysym(&report.xkey, 0) == XK_q){ // quit           
                    finished=true;
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_d){ // increment depth    
                	depth++;
                	draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_e){ // decrement depth          
					depth--;
					draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Left){          
                    z=z-step;
                    az = abs(z);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Right){           
                    z=z+step;
                    az = abs(z);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Up){          
                    z=z+(step*I);
                    az = abs(z);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Down){           
                    z=z-(step*I);
                    az = abs(z);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_j){          
                    w=w+step;
                    aw = abs(w);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_l){           
                    w=w-step;
                    aw = abs(w);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_i){          
                    w=w-(step*I);
                    aw = abs(w);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_k){           
                    w=w+(step*I);
                    aw = abs(w);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_s){ // speed up       
                    step=step*10.0;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_a){ // slow down         
                    step=step/10.0;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_n){ // increase mesh       
                    mesh=mesh*2;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_m){ // decrease mesh       
                    if(mesh>1){
	                    mesh=mesh/2;
    	                draw();
    	            };
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_g){ // zoom out       
                    wind=wind*2.0;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_b){ // mode toggle 
                	  s=sync;	// remember sync mode
                	  reinitialize(center,center);
                	  sync=s;	// don't reset sync mode
                    mode=1-mode;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_v){ // sync cycle 
                    sync=(sync+1)%3;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_c){ // toggle color_ifs
                    color_ifs=1-color_ifs;
                    draw();
                    break;
                } else if (XLookupKeysym(&report.xkey, 0) == XK_f) { //toggle disconnected depth
                    disconnection_depth = !disconnection_depth;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_t){ 
                    draw_trap_mode = !draw_trap_mode;
                    draw();
                    break;
		            } else if (XLookupKeysym(&report.xkey, 0) == XK_x) { //toggle chunky ifs
                    chunky_ifs = !chunky_ifs;
                    if (mode == 0) draw();
                    break;
                }

            default:
            	break;
        }
    }
}
