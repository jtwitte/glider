#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "distance.h"
#include "utils.h"
#include "lambda2.h"
//#include "tracer.h"

scalar f[], d[]; // make a global field
//scalar trac[];
vertex scalar phi[];
face vector s[];

int    LEVEL   = 12; 		// define grid resolution
double uemax   = 0.02; 		// define mesh adaption criterion
double MU0 = 1.08e-3;//3.258e-3; 		// dynamic viscosity
/*double maxruntime = 1438;  just less than 24 hours in minutes */
double maxruntime = 700; 	// just less than 12 hours in minutes
double WIDTH = 5;
double MONIT = 2.45;

//scalar * tracers = {trac};

/* boundary conditions for velocity with unit velocity on the left-hand side and outflow on the right-hand side  */

u.n[top]  = dirichlet(0.33);
p[top]    = neumann(0.);
pf[top]   = neumann(0.);
//trac[left]  = dirichlet(y<0);
//trac[left]  = dirichlet(0.33);

u.n[bottom] = neumann(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);

int main (/*int argc, char *argv[]*/) {

/* kinametic viscosity of water 1.05x10^-6 m^2 s^-1 at 20 degrees*/
/* kinametic viscosity of water 1.83x10^-6 m^2 s^-1 at 0 degrees*/
/* dynamic viscosity (MU0) of water 1.08 × 10−3 Pa s at 20 degrees*/
/* dynamic viscosity (MU0) of water 1.88 × 10−3 Pa s at 0 degrees*/
/* http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html */
/* actual Re = l*V/v , 1.8*0.33./1.05x10^-6 = 5.65714e5 */
/* use v= 1.8*1./ 5.65714e5 = 3.1818*10^-6 */
/* use MU0 = v * density = 3.1818*10^-6 * 1024 = 3.258e-3 */
  const face vector muc[] = { 1.05e-6 , 1.05e-6, 1.05e-6 }; // { 3.1818e-6, 3.1818e-6, 1.1818e-6 };
  mu = muc;

  run();
}

/* setup of geometry and domain*/

event init (t = 0) {
  if (!restore (file = "restart")) {
  fprintf(ferr,"initialising\n");

/*if (!restore (file = "restart")) {
    refine (sq(x) + sq(y - Z0) + sq(z) - sq(0.75) < 0 && level < LEVEL);
    fraction (f, sq(x) + sq(y - Z0) + sq(z) - sq(.5));
  }
*/
  coord * p = input_stl (fopen ("sphere.stl","r"));//"omg_glider_final_2_5degree_trans_flip_binary.stl", "r"));
/* coord min, max;
   bounding_box (p, &min, &max);  
   double maxl = -HUGE;
   foreach_dimension()
     if (max.x - min.x > maxl)
       maxl = max.x - min.x;
*/  

  init_grid (8);
  size (WIDTH);
/*  fprintf(ferr,"%f  (max.x+min.x)/2. - L0/2. \n", (max.x+min.x)/2. - L0/2.);
    fprintf(ferr,"%f  (max.x+min.y)/2. - L0/2. \n", (max.y+min.y)/2. - L0/2);
    fprintf(ferr,"%f  (max.z+min.z)/2. - L0/2. \n", (max.z+min.z)/2. - L0/2);
 
   origin ((max.x + min.x)/2. - L0/2,          
 * 	    (max.y + min.y)/2. - L0/2,
 *  	    (max.z + min.z)/2. - L0/2); */
  origin( -L0/2, -L0/2., -L0/2. );
  distance (d, p);

  while (adapt_wavelet ({d}, (double[]){1e-4}, LEVEL).nf);
 
  foreach()
  u.x[] = 0;//.33; // 1.;
                }
  else {
  fprintf(ferr,"restarting\n"); 
  }
}
 
/* calculating ? */
 
event glider (i++) {
  fprintf(ferr,"define zero velocity inside geometry\n"); 
  coord vc = {0.,0.,0.};			// the velocity of the glider
 
/* recalculate the fraction field at each time step or the geometry goes gets very noisy*/

  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f, s);
  
  foreach() 
    foreach_dimension() 
        u.x[] = f[]*vc.x + (1. - f[])*u.x[];
  boundary ((scalar *){f,u});			// include f as a boundary ?
 
}


event snapshot (t = 0; t += 1; t <= 25) {
 
  fprintf(ferr,"saving file\n");
  dump ( file = "restart" );

  char name[80];
  sprintf (name, "snapshot-%g-%d.gfs", t,LEVEL);
 
  scalar pid[];
  scalar l2[], vyz[] ;
  scalar omega[];
 
  foreach() { 
    vyz[] = ( (u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
    pid[] = fmod(pid()*(npe() + 37), npe()); };
 
  lambda2 (u, l2);
  vorticity (u, omega);

  boundary ({pid,vyz});
  
  output_gfs (file = name ); 
 
}

/* adaptive grid size*/

event adapt (i++) {
	astats st= adapt_wavelet ((scalar *){u}, (double[]){uemax,uemax,uemax},LEVEL);
  	fprintf (ferr, "  flow field # grid %ld cells, refined %d cells, coarsened %d cells\n", grid->tn, st.nf, st.nc);
}

/* print logfile */
 
event logfile (i++){
	fprintf(stderr,"i= %d t= %.4f dt=%g mg[ %d %d %d ] max(u,v,w)= %.3f %.3f %.3f perf(t,speed)= %.3f %.3f \n",i, t, dt, mgp.i, mgpf.i, mgu.i, statsf(u.x).max, statsf(u.y).max, statsf(u.z).max, perf.t, perf.speed);     

	FILE * timefile = fopen ("scaling.dat", "aw");
	fprintf (timefile, "%d %.4f %.3f %.3f \n", i, t, perf.t, perf.speed);
	fclose(timefile);
}


event timeseries(i++) {

// MONITORING LOCATIONS
    monitor(i,t,{p,u}, -MONIT, MONIT, 0,"inlet_L5.dat");
    monitor(i,t,{p,u}, 0, MONIT, 0,"opposite_L5.dat");
    monitor(i,t,{p,u}, MONIT, -MONIT, 0,"outlet_L5.dat");
    monitor(i,t,{p,u}, 0, 0, MONIT,"above_L5.dat");
    monitor(i,t,{p,u}, MONIT, 0, 0,"inlet_front_L5.dat");
     
//added in 0.005 (5mm in x) to probe monitoring points so they are not on the surface of the geometry
    monitor(i,t,{p,u}, -0.79281, 0.02224, 0.18038,"probe1.dat");
    monitor(i,t,{p,u}, -0.79268, 0.02199, 0.20597,"probe2.dat");
    monitor(i,t,{p,u}, -0.78936, 0.00000, 0.21854,"probe3.dat");
    monitor(i,t,{p,u}, -0.75805, -0.02199, 0.20447,"probe4.dat");
    monitor(i,t,{p,u}, -0.74699, 0.00000, 0.16585,"probe5.dat");
           
//added in 0.005 (5mm in z) to Oxygen to lift of surface
    monitor(i,t,{p,u}, 1.0851, -0.03200, -0.00768,"oxygen.dat");
    monitor(i,t,{p,u}, -0.03914, -0.15324, -0.05431,"ctd.dat");
              
    monitor(i,t,{p,u}, 1.15, 0, -0.0393,"wake1_centre.dat");
    monitor(i,t,{p,u}, 1.35, 0, -0.0393,"wake2_centre.dat");
    monitor(i,t,{p,u}, 1.55, 0, -0.0393,"wake3_centre.dat");
    monitor(i,t,{p,u}, 1.75, 0, -0.0393,"wake4_centre.dat");
    monitor(i,t,{p,u}, 1.95, 0, -0.0393,"wake5_centre.dat");
    monitor(i,t,{p,u}, 2.15, 0, -0.0393,"wake6_centre.dat");
    monitor(i,t,{p,u}, 2.35, 0, -0.0393,"wake7_centre.dat");
                   
    monitor(i,t,{p,u}, 1.15, 0, 0.0166,"wake1_top.dat");
    monitor(i,t,{p,u}, 1.35, 0, 0.0166,"wake2_top.dat");
    monitor(i,t,{p,u}, 1.55, 0, 0.0166,"wake3_top.dat");
    monitor(i,t,{p,u}, 1.75, 0, 0.0166,"wake4_top.dat");
    monitor(i,t,{p,u}, 1.95, 0, 0.0166,"wake5_top.dat");
    monitor(i,t,{p,u}, 2.15, 0, 0.0166,"wake6_top.dat");
    monitor(i,t,{p,u}, 2.35, 0, 0.0166,"wake7_top.dat");
                                
    monitor(i,t,{p,u}, 1.15, 0, -0.0952,"wake1_lower.dat");
    monitor(i,t,{p,u}, 1.35, 0, -0.0952,"wake2_lower.dat");
    monitor(i,t,{p,u}, 1.55, 0, -0.0952,"wake3_lower.dat");
    monitor(i,t,{p,u}, 1.75, 0, -0.0952,"wake4_lower.dat");
    monitor(i,t,{p,u}, 1.95, 0, -0.0952,"wake5_lower.dat");
    monitor(i,t,{p,u}, 2.15, 0, -0.0952,"wake6_lower.dat");
    monitor(i,t,{p,u}, 2.35, 0, -0.0952,"wake7_lower.dat");
                                    
    monitor(i,t,{p,u}, 1.15, 0, 0.150,"wake1_upper.dat");
    monitor(i,t,{p,u}, 1.35, 0, 0.150,"wake2_upper.dat");
    monitor(i,t,{p,u}, 1.55, 0, 0.150,"wake3_upper.dat");
    monitor(i,t,{p,u}, 1.75, 0, 0.150,"wake4_upper.dat");
    monitor(i,t,{p,u}, 1.95, 0, 0.150,"wake5_upper.dat");
    monitor(i,t,{p,u}, 2.15, 0, 0.150,"wake6_upper.dat");
    monitor(i,t,{p,u}, 2.35, 0, 0.150,"wake7_upper.dat");

} 


event force (i++) {
  
  vector dux[]; vector duy[]; vector duz[];
  vector Sx[]; vector Sy[]; vector Sz[];

foreach(){
    //dux.x[] = (u.x[1,0] - u.x[-1,0])/2./Delta;
    dux.x[] = (u.x[1,0,0] - u.x[-1,0,0])/2./Delta;
    dux.y[] = (u.x[0,1,0] - u.x[0,-1,0])/2./Delta;
    dux.z[] = (u.x[0,0,1] - u.x[0,0,-1])/2./Delta;

    duy.y[] = (u.y[0,1,0] - u.y[0,-1,0])/2./Delta;
    duy.x[] = (u.y[1,0,0] - u.y[-1,0,0])/2./Delta;
    duy.z[] = (u.y[0,0,1] - u.y[0,0,-1])/2./Delta;

    duz.z[] = (u.z[0,0,1] - u.z[0,0,-1])/2./Delta;
    duz.x[] = (u.z[1,0,0] - u.z[-1,0,0])/2./Delta;
    duz.y[] = (u.z[0,1,0] - u.z[0,-1,0])/2./Delta;
}    

foreach(){
    Sx.x[] = -p[] + MU0*(2*dux.x[]);
    Sx.y[] = MU0*(dux.y[] + duy.x[]);
    Sx.z[] = MU0*(dux.z[] + duz.x[]);

    Sy.y[] = -p[] + MU0*(2*duy.y[]);
    Sy.x[] = MU0*(duy.x[] + dux.y[]);
    Sy.z[] = MU0*(duy.z[] + duz.y[]);

    Sz.z[] = -p[] + MU0*(2*duz.z[]);
    Sz.x[] = MU0*(duz.x[] + dux.z[]);
    Sz.y[] = MU0*(duz.y[] + duy.z[]);

}

  double Faerox=0; double Faeroy=0; double Faeroz=0;
  double Strainxx=0; double Strainyx=0; double Strainzx=0;
  double Strainxy=0; double Strainyy=0; double Strainzy=0;
  double Strainxz=0; double Strainyz=0; double Strainzz=0;
  /*double Faero=0;*/
  foreach(reduction(+:Faerox) reduction(+:Faeroy) reduction(+:Faeroz)){

    if (f[] > 1e-4 && f[] < 1. - 1e-4) {
    //if (f[] < 1. + 1e-4 && f[] > 1. - 1e-4) {
        coord m = mycs (point, f);
        double alpha = plane_alpha (f[], m);
        coord pp;

        Faerox += m.x*(Sx.x[]+Sy.x[]+Sz.x[])*pow(Delta, dimension-1)*plane_area_center(m, alpha, &pp);
        Faeroy += m.y*(Sx.y[]+Sy.y[]+Sz.y[])*pow(Delta, dimension-1)*plane_area_center(m, alpha, &pp);
        Faeroz += m.z*(Sx.z[]+Sy.z[]+Sz.z[])*pow(Delta, dimension-1)*plane_area_center(m, alpha, &pp);
	
	Strainxx += Sx.x[];
	Strainyx += Sy.x[];
	Strainzx += Sz.x[];
	Strainxy += Sx.y[];
	Strainyy += Sy.y[];
	Strainzy += Sz.y[];
	Strainxz += Sx.z[];
	Strainyz += Sy.z[];
	Strainzz += Sz.z[];

    }
}

  //Faero = sqrt(sq(Faerox) + sq(Faeroy) + sq(Faeroz));

  //fprintf (ferr, " %f %d  %f %f %f\n", t, i, Faerox, Faeroy, Faeroz);
  //float forces[3];
  float forces[12];
  FILE * ldrag = fopen ("forces.dat", "aw");
  forces[0] = Faerox;
  forces[1] = Faeroy;
  forces[2] = Faeroz;
 
  forces[3] = Strainxx;
  forces[4] = Strainyx;
  forces[5] = Strainzx;
  forces[6] = Strainxy;
  forces[7] = Strainyy;
  forces[8] = Strainzy;
  forces[9] = Strainxz;
  forces[10] = Strainyz;
  forces[11] = Strainzz;
 
  if (pid() == 0){ // master
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, &forces, 12, MPI_DOUBLE, MPI_MIN, 0,
                    MPI_COMM_WORLD);
 
      @endif
        fprintf (ferr, " %f %d  %f %f %f %f %f %f %f %f %f %f %f %f \n", t, i, forces[0], forces[1], forces[2], forces[3], forces[4], forces[5], forces[6], forces[7], forces[8], forces[9], forces[10], forces[11]);
        fprintf (ldrag, " %f %d  %f %f %f %f %f %f %f %f %f %f %f %f \n", t, i, forces[0], forces[1], forces[2], forces[3], forces[4], forces[5], forces[6], forces[7], forces[8], forces[9], forces[10], forces[11]);
    }
  @if _MPI
  else // slave
    MPI_Reduce (&forces, NULL, 12, MPI_DOUBLE, MPI_MIN, 0,
                  MPI_COMM_WORLD);
  @endif
 
  fclose(ldrag);
	
}
 
/* not working with mpi yet
event timeseries(i++) {
   fprintf(ferr,"monitoring point\n");
    FILE * fp = fopen ("test.dat", "w");
     fprintf (fp, " time " " iteration " " pressure " " u comp" " v comp " " w comp\n" );
     float p_point = interpolate (p  , 1.4, 1 , 1);
     float u_point = interpolate (u.x, 1.4, 1 , 1);
     float v_point = interpolate (u.y, 1.4, 1 , 1);
     float w_point = interpolate (u.z, 1.4, 1 , 1);
     fprintf (fp, " %f %d  %f %f %f %f \n", t, i, p_point, u_point,v_point, w_point );
     fprintf (ferr, " %f %d  %f %f %f %f \n", t, i, p_point, u_point,v_point, w_point );
}*/

event runtime (i += 1) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t/60 >= maxruntime) {
    dump (file = "restart"); 			// so that we can restart
    return 1; 					// exit
  }
}


