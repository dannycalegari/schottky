/* interface.cc */

void ifs::user_interface(){
	bool finished;
	point p,pp;
	cpx mz;
	stringstream T;
	int s;
	
	finished=false;
	
	cout << "Entering interactive mode. Press [q] to exit.\n";
	cout << "(note: mouse must be focussed in X-window)\n";
	while(finished==false){
		XNextEvent(display, &report);
		switch(report.type) {
			case ButtonPress:
				p=mouse_location();
				if(mode==1){	// mandelbrot mode
					zoom(p);
					draw();
				};
				break;
			case MotionNotify:
				p=mouse_location();
				if(mode==1){	// mandelbrot mode
					mz=point_to_cpx(p);
					mz=(mz*wind)+center;
					p.x=1200;	// 1200
					p.y=500;
					T.str("");
					T << "mouse location " << mz;
					pp=p;
					pp.y=pp.y+100;
					draw_box(pp, 400, 0xFFFFFF);
					draw_text(p, T, 0x000000);
				};
				break;
			case KeyPress:
				if(XLookupKeysym(&report.xkey, 0) == XK_q){ // quit           
                    finished=true;
                    XCloseDisplay(display);
                    exit(0);
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
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Right){           
                    z=z+step;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Up){          
                    z=z+(step*I);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_Down){           
                    z=z-(step*I);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_j){          
                    w=w+step;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_l){           
                    w=w-step;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_i){          
                    w=w-(step*I);
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_k){           
                    w=w+(step*I);
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
                	initialize(center,center);
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
                };

            default:
            	break;
        };
    };
};
