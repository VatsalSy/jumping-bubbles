/* Title: Saving images with bview
# Authors: Vatsal & Youssef
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];
char filename[80];

trace
void output_facets_v2 (scalar c, FILE * fp = stdout, face vector s = {{-1}})
{
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (c[], n);
#if dimension == 1
      fprintf (fp, "%g\n", x + Delta*alpha/n.x);
#elif dimension == 2
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
	fprintf (fp, "%g %g\n%g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets (n, alpha, v, 1.1);
      for (int i = 0; i < m; i++)
	fprintf (fp, "%g %g %g\n",
		 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	fputc ('\n', fp);
#endif
    }

  fflush (fp);
}

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);

  f[bottom] = dirichlet(1.);
  restore (file=filename);

  #if TREE
  // make sure we prolongate properly
  void (* prolongation) (Point, scalar) = f.prolongation;
  if (prolongation != fraction_refine) {
    f.prolongation = fraction_refine;
    f.dirty = true;
  }
  #endif // TREE


  output_facets_v2(f, ferr);

}
