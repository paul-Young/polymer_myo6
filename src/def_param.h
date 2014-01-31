#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structure.h"

#define tot_amino  13     /* number of amino acids -- needed to insure that structures 1 and 2 are compatible */ 
#define tot_pol1   7   

#define h        0.04     /* integration step */
#define zeta     25.0     /* (zeta*(1/tau_L)*(eh/temp) sets the natural time for the overdamped motion */
#define a0       3.8 	    /* soft-core repulsive distance */
#define kspring_cov 200.0  /* spring constant in the harmonic bond between adjacent beads */
#define el       20.0     /* energy scale for purely repulsive LJ pot. */
#define k_bend	 14.0      /* bending potential strength - define persistence length */
#define R_limit  10.0     /* max distance between beads before program exits, compare to a0 */
#define stepsim  10000000000 /* number of simulation steps for heating under zero applied force */
#define nav      10000    /* frequency to get information about chain conformation in update */
#define nav1     nav      /* frequency to save data in out file  */
#define nav2     (100*nav)   /* frequency to save chain conformation in Struct_data */
#define temp     0.6      /* temperature (in units of (eh/kB)); 0.6 kcal/mol = 300 K */
#define eps      1.E-6
#define INF      100000.0
#define R_sigma  10.0      /* ratio of bead distance to soft core repulsion distance, above this, the repulsion is ignored*/
//#define vc 	 1.0 /* end-tangent constraint strength - can play */
//#define thetac	 30.0 /* end-tangent constraint angle with respect to actin or the z-axis - can play */
#define DR2_limit 1.0 /* capture radius^2 */
#define xb	0.0
#define yb	0.0
#define zb	36.0	// capture site 
#define xbb	0.0
#define ybb	0.0
#define zbb	-36.0	// backward binding site
//#define pull_force	1.0  // pulling force on the dimerization site

double vc,thetac,pull_force;
double dr2_min; // min distance to binding site
int bound; // flag for reaching binding site
double tangx, tangy, tangz; // end tangent constaint related

long long int step;      // to acommodate long runs
int run; // for multiple runs (must be an interger less than 10000, because run-10000 is used as random seed)
int    *order_N, *nr_N;
double **R0;
double epot_bend,epot_bond,epot_rep,epot_tangent,tempav,epot,rg;
int    mseed;

atom *Atom;
atom *Amino, *d_Amino, *Amino0;

FILE *fout, *fout1, *fcd, *fdata, *f_ext;
char filename_b[100];
char   *str_mseed;
char   **type_N, **mol_N;

struct coord{
  double x, y, z;
};
struct coord c;
