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
                } else if (XLookupKeysym(&report.xkey, 0) == XK_p) { //toggle drawing 1/2
                    draw_contains_half = !draw_contains_half;
                    draw();
                    break;
                } else if(XLookupKeysym(&report.xkey, 0) == XK_w){ 
                    find_close_uv_words = !find_close_uv_words;
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
                } else if (XLookupKeysym(&report.xkey, 0) == XK_r) { //toggle trap-like vectors
                    find_trap_like_vectors = !find_trap_like_vectors;
                    draw();
                    break;
                } else if (XLookupKeysym(&report.xkey, 0) == XK_o) { //start drawing a loop
                    if (mode == 1) { //only do it in mandlebrot mode
                      //get a loop
                      std::vector<cpx> loop(0);
                      input_loop(loop);
                      find_TLB_along_loop(loop, true, 1);
                      break;
                    }
                }

            default:
            	break;
        }
    }
}



//in mandlebrot mode, read in a sequence of points
void ifs::input_loop(std::vector<cpx>& loop) {
  bool done_drawing = false;
  Point2d<int> p;
  int wcol = X.get_rgb_color(1,1,1);
  std::stringstream T;
  Point2d<int> pp;
  
  p.x = drawing_width + 15;
  p.y = drawing_width/2 + 50;
  T.str("");
  T << "*click to draw path/loop*";
  X.draw_text(p, T, 0);
  p.y-=20;
  T.str("");
  T << "[o] to finish";
  X.draw_text(p, T, 0);
  
  while (!done_drawing) {
	  XEvent report;
	  X.get_next_event(report);
	  if (report.type == MotionNotify) {
	    p = X.mouse_location();
      cpx mz = point_to_cpx(p);
      mz=(mz*wind)+center;
      p.x=drawing_width + 15;
      p.y=drawing_width/2;
      T.str("");
      T << "mouse location " << mz;
      pp=p;
      pp.y=pp.y-30;
      X.draw_filled_rectangle(pp, 300, 40, 0xFFFFFF);
      X.draw_text(p, T, 0x000000);
    
    } else if (report.type == ButtonPress) {
      p = X.mouse_location();
      cpx mz = point_to_cpx(p);
      mz = (mz*wind)+center;   //now mz is the real complex number
      X.draw_dot(p, wcol);
      if (loop.size() > 0) { 
        X.draw_line(cpx_to_point_mandlebrot(loop.back()), p, wcol);
      }
      loop.push_back(mz);
    
    } else if (report.type == KeyPress && 
               XLookupKeysym(&report.xkey, 0) == XK_o) {
      X.draw_line(cpx_to_point_mandlebrot(loop.back()), 
                  cpx_to_point_mandlebrot(loop.front()), wcol);
      done_drawing= true;
    }
  }
  
  p.x = drawing_width + 15;
  p.y = drawing_width/2 + 20;
  X.draw_filled_rectangle(p, 300, 50, wcol);
	    
	
}










































