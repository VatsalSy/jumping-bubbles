/* Title: Interpolating data from dump files: gfs2oogl style
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

char filename[80];
int ny, nz;
double ymin, zmin, ymax, zmax, xSlice, Oh;
bool linear;
scalar * list = NULL;

scalar vel[], D2c[];

#define Mu21 (1.00e-3)
#define Rho21 (1.00e-3)

int main(int a, char const *arguments[]){

  // boundary conditions
  u.t[bottom] = dirichlet(0.);
  u.r[bottom] = dirichlet(0.);
  f[bottom] = dirichlet(1.);

  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary ({f,u.x,u.y,u.z});

  ymin = HUGE;
  foreach_boundary(left){
    if (y < ymin) ymin = y;
  }

  ymax = atof(arguments[2]);

  zmin = 0.0;

  zmax = atof(arguments[3]);

  xSlice = atof(arguments[4]);

  nz = atoi(arguments[5]);

  Oh = atof(arguments[6]);

  linear = arguments[7];

  rho1 = 1.0; mu1 = Oh;
  rho2 = Rho21; mu2 = Mu21*Oh;

  list = list_add (list, f);
  // list = list_add (list, u.y);
  // list = list_add (list, u.z);
  list = list_add (list, vel);
  list = list_add (list, D2c);

  double DeltaMin = HUGE;

  foreach(){
    vel[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));

    double D2 = 0.;
    foreach_dimension(){
      double DII = (u.x[1,0,0]-u.x[-1,0,0])/(2*Delta);
      double DIJ = 0.5*((u.x[0,1,0]-u.x[0,-1,0] + u.y[1,0,0] - u.y[-1,0,0])/(2*Delta));
      double DIK = 0.5*((u.x[0,0,1]-u.x[0,0,-1] + u.z[1,0,0] - u.z[-1,0,0])/(2*Delta));
      D2 += sq(DII) + sq(DIJ) + sq(DIK);
    }
    D2c[] = 2*(mu(f[]))*D2; //*cube(Delta);
    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10.);
    } else {
      D2c[] = -10.;
    }

    DeltaMin = DeltaMin > Delta ? Delta : DeltaMin;

  }

  if ((linear == true) && (DeltaMin < 4*((double)((zmax-zmin)/(nz))))){
    linear = false;
  }

  if (linear == false){
    FILE * fp = ferr;
    foreach_boundary(left){
      // if (y-2*Delta > ymin && y+2*Delta < ymax && z-2*Delta > zmin && z+2*Delta < zmax)
        fprintf (fp, "%g %g %g %g %g\n", y, z, f[], vel[], D2c[]);
    }
  } else {
    FILE * fp = ferr;
  
    double DeltaZ = (double)((zmax-zmin)/(nz));
    int ny = (int)((ymax-ymin)/DeltaZ);
    double DetlaY = (double)((ymax-ymin)/(ny));

    // fprintf (ferr, "ny = %d, nz = %d\n", ny, nz);
    
    int len = list_len(list);
    double ** field = (double **) matrix_new (ny, nz, len*sizeof(double));
    for (int i = 0; i < ny; i++) {
      double y = DetlaY*i + ymin;
      for (int j = 0; j < nz; j++) {
        double z = DeltaZ*j + zmin;
        int k = 0;
        for (scalar s in list){
          field[i][len*j + k++] = interpolate (s, xSlice, y, z);
        }
      }
    }

    for (int i = 0; i < ny; i++) {
      double y = DetlaY*i + ymin;
      for (int j = 0; j < nz; j++) {
        double z = DeltaZ*j + zmin;
        fprintf (fp, "%g %g", y, z);
        int k = 0;
        for (scalar s in list){
          fprintf (fp, " %g", field[i][len*j + k++]);
        }
        fputc ('\n', fp);
      }
    }
    fflush (fp);
    matrix_free (field);
  }
}