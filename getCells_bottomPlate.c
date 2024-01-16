/* Title: Saving images with bview
# Authors: Vatsal & Youssef
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  restore (file=filename);
  foreach_boundary(bottom){
    fprintf (ferr, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
	       x - Delta/2., y, z - Delta/2.,
	       x - Delta/2., y, z + Delta/2.,
	       x + Delta/2., y, z + Delta/2.,
	       x + Delta/2., y, z - Delta/2.);
  }
}
