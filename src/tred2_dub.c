#include <stdio.h>
#include <math.h>

void tred2(double a[4][4], int n, double *d, double *e)
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  /*for(i = 1 ; i <= n ; i ++){
    for(j = 1 ; j <= n ; j ++){
      printf("in tred2: %3d %3d : %10.5f\n", i, j, a[i][j]);
    }
    }*/

  for (i = n ; i >= 2 ; i --) {
    l = i - 1;
    h = scale = 0.0;

    if (l > 1) {

      for (k = 1 ; k <= l ; k ++){

	scale += fabs(a[i][k]);

	//printf("%5d %5d a = %g scale = %g\n", i, k, a[i][k], scale);
      
      }

      //printf("i = %5d l = %5d: scale = %g\n", i, l, scale);

      if (scale == 0.0){
       	e[i] = a[i][l];

	//printf("\n\nFor i = %5d scale = %g\n", i, scale);

      }

      else {

	for (k = 1 ; k <= l ; k ++) {

	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	
	  //printf("%5d %5d a = %g h = %g\n", i, k, a[i][k], h);
	
	}

	f = a[i][l];
	g = f>0 ? -sqrt(h) : sqrt(h);
	e[i] = scale * g;
	h -= f * g;
	a[i][l] = f - g;
	f = 0.0;

	for (j = 1 ; j <= l ; j ++) {
	  /* Next statement can be omitted if eigenvectors not wanted */
	  a[j][i] = a[i][j]/h;
	  g = 0.0;

	  for(k = 1 ; k <= j ; k ++)
	    g += a[j][k] * a[i][k];

	  for (k = (j+1) ; k <= l ; k ++)
	    g += a[k][j]*a[i][k];

	  e[j] = g/h;
	  f += e[j]*a[i][j];
	}

	hh = f/(h+h);

	for(j = 1 ; j <= l ; j ++) {

	  f = a[i][j];
	  e[j] = g = e[j] - hh * f;

	  for (k = 1 ; k <= j ; k ++){
	    a[j][k] -= (f * e[k] + g * a[i][k]);
	  
	    //printf("%5d %5d a = %g\n", j, k, a[j][k]);

	  }

	}

      }/* end scale is non-zero */

    } /* end l > 1 */

    else
      e[i] = a[i][l];

    d[i] = h;

  }
  /* end cycle "i" */


  /* Next statement can be omitted if eigenvectors not wanted */
  d[1]=0.0;
  e[1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i = 1 ; i <= n ; i ++) {
    l = i-1;
    if (d[i]) {
      for (j = 1 ; j <= l ; j ++) {
	g = 0.0;
	for (k = 1 ; k <= l ; k ++)
	  g += a[i][k]*a[k][j];
	for (k = 1 ; k <= l ; k ++)
	  a[k][j] -= g*a[k][i];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for (j = 1 ; j <= l ; j ++) a[j][i]=a[i][j]=0.0;
  }

}

