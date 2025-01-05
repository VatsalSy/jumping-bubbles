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
int ny, nx;
double ymin, xmin, ymax, xmax, zSlice, Oh;
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
  foreach_boundary(back){
    if (y < ymin) ymin = y;
  }

  ymax = atof(arguments[2]);

  xmin = 0.0;

  xmax = atof(arguments[3]);

  zSlice = atof(arguments[4]);

  nx = atoi(arguments[5]);

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

  if ((linear == true) && (DeltaMin < 4*((double)((xmax-xmin)/(nx))))){
    linear = false;
  }

  if (linear == false){
    FILE * fp = ferr;
    foreach_boundary(back){
      // if (y-2*Delta > ymin && y+2*Delta < ymax && x-2*Delta > xmin && x+2*Delta < xmax)
        fprintf (fp, "%g %g %g %g %g\n", y, x, f[], vel[], D2c[]);
    }
  } else {
    
    FILE * fp = ferr;

    double DeltaX = (double)((xmax-xmin)/(nx));
    int ny = (int)((ymax-ymin)/DeltaX);
    double DetlaY = (double)((ymax-ymin)/(ny));

    // fprintf (ferr, "ny = %d, nz = %d\n", ny, nz);
    
    int len = list_len(list);
    double ** field = (double **) matrix_new (ny, nx, len*sizeof(double));
    for (int i = 0; i < ny; i++) {
      double y = DetlaY*i + ymin;
      for (int j = 0; j < nx; j++) {
        double x = DeltaX*j + xmin;
        int k = 0;
        for (scalar s in list){
          field[i][len*j + k++] = interpolate (s, x, y, zSlice);
        }
      }
    }

    for (int i = 0; i < ny; i++) {
      double y = DetlaY*i + ymin;
      for (int j = 0; j < nx; j++) {
        double x = DeltaX*j + xmin;
        fprintf (fp, "%g %g", y, x);
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