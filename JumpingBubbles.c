/* Title: Jumping Drops
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "distance.h"

#define MINlevel 2                                              // maximum level

#define tsnap (1e-2)
#define tsnap2 (1e-4)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity

#define Mu21 (1.00e-3)
#define Rho21 (1.00e-3)

// domain
#define Ldomain 4                                // Dimension of the domain

// boundary conditions
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
f[bottom] = dirichlet(1.);

double tmax, Oh;
int MAXlevel;                                              // maximum level

int main() {

  tmax = 1e-2;
  Oh = 0.01; // <0.001/sqrt(1000*0.072*0.001)>
  MAXlevel = 9;

  init_grid (1 << MINlevel);
  L0=Ldomain;
  fprintf(ferr, "tmax = %g. Oh = %g\n",tmax, Oh);
  rho1 = 1.0; mu1 = Oh;
  rho2 = Rho21; mu2 = Mu21*Oh;
  f.sigma = 1.0;

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    char filename[60];
    sprintf(filename,"InitialCondition.stl");
    FILE * fp = fopen (filename, "r");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord * p = input_stl (fp);
    fclose (fp);
    coord min, max;

    bounding_box (p, &min, &max);
    fprintf(ferr, "xmin %g xmax %g\nymin %g ymax %g\nzmin %g zmax %g\n", min.x, max.x, min.y, max.y, min.z, max.z);
    fprintf(ferr, "x0 = %g, y0 = %g, z0 = %g\n", 0., - 1 - L0/pow(2,MAXlevel+1), (min.z+max.z)/2.);
    origin (0., - 1 - L0/pow(2,MAXlevel), (min.z+max.z)/2.);

    scalar d[];
    distance (d, p);
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-6, 1e-6*L0}, MAXlevel).nf);
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
  	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    }
    fractions (phi, f);

    foreach () {
      foreach_dimension(){
        u.x[] = 0.0;
      }
    }

    dump (file = "dumpInit");
    fprintf(ferr, "Done with initial condition!\n");
    // return 1;
  }
}

event adapt(i++) {
  scalar KAPPA[];
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, KAPPA, u.x, u.y, u.z},
     (double[]){fErr, KErr, VelErr, VelErr, VelErr},
      MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 0.5*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]))*clamp(1.-f[], 0., 1.)*cube(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

}
