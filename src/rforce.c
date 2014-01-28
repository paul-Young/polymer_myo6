
/* this is the program "rforce.c". It calculates the random force in Langevin simulations. */

#include "def_param.h"

double gasdev(), ran2();

/* the random force is white noise (i.e., it is Gaussian distributed) */
void rforce(){

  int i;
  double var;

  var = sqrt(2.0 * temp/h);

  for(i = 1 ; i <= tot_amino ; i ++){
    Amino[i].fcx = var * gasdev(&mseed) * sqrt(Amino0[i].zeta1);
    Amino[i].fcy = var * gasdev(&mseed) * sqrt(Amino0[i].zeta1);
    Amino[i].fcz = var * gasdev(&mseed) * sqrt(Amino0[i].zeta1);
//    printf("random forces: %f, %f, %f \n", Amino[i].fcx, Amino[i].fcy, Amino[i].fcz);
  }
}


double gasdev(int *idum){

  static int   iset = 0;
  static double gset;
  double        fac, rsq, v1, v2;
  
  if(iset == 0){
    do{
      v1 = 2.0 * ran2(idum) - 1.0;
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1*v1+v2*v2;

    } while( rsq >= 1.0 || rsq == 0.0);
    
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;

    return v2*fac;
  }

  else{

    iset = 0;

    return gset;
  }

}

                 /* this is the program ran2() */

/* it generates a random double number b/w 0 and 1. It is set-up for a DEC alpha -- 64-bit -- (that is, 
idum2, iy and iv are static integers; for a pentium -- 32-bit -- we need to make them long integers). */


#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1 -1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1 + IMM1/NTAB)
#define EPS1   1.2e-7
#define RNMX  (1.0 - EPS1) 

double ran2(int *idum){
  int         j;
  int         k;
  static int  idum2 = 123456789;
  static int  iy = 0;
  static int  iv[NTAB];
//  long int  idum2 = 123456789;
//  long int  iy = 0;
//  long int  iv[NTAB];
  double      tempran;

  if(*idum <= 0){   /* initialize */
  
    if(-(*idum) < 1) *idum = 1; /* be sure to prevent idum = 0 */
    else *idum = - (*idum);
    idum2 = (*idum);

    for(j = NTAB+7 ; j >= 0 ; j --){ /*load shuffle table (after 8 warm-ups) */
      k = (*idum)/IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum += IM1;
      if(j < NTAB) iv[j] = *idum;
    }

    iy = iv[0];

  }

  k = (*idum)/IQ1;                    /* start here when not init. */
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((tempran = AM*iy) > RNMX) return RNMX;
  else return tempran;

}

