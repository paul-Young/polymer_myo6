/* this is the program initialize.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "def_param.h"


int initialize(){

  int i,j;
  double dx, dy, dz, dr;

  printf("....initializing...structure 1 \n"); 
  R0      = malloc((tot_amino+1) * sizeof(double *)); 
  nr_N    = (int*)malloc((tot_amino+1) * sizeof(int));
  mol_N   = (char **)malloc((tot_amino+1) * sizeof(char*));

  for(i = 1 ; i <= tot_amino ; i ++){
    mol_N[i]  = (char *)malloc((tot_amino+1) * sizeof(char));
  }
  
  
  for(i = 1 ; i <= tot_amino ; i ++){ 
    R0[i]     = malloc((tot_amino+1) * sizeof(double));
  }

  /* build R0 matrices for Fene potential */
  for(i = 1 ; i <= (tot_amino-2) ; i ++){
    for(j = (i+2) ; j <= tot_amino ; j ++){
      	R0[i][j]= Amino[i].a1 + Amino[j].a1;    // set a0 to be the repulsive radius 
    }
  }
  for(i = 1 ; i <= (tot_amino-1) ; i ++){    
    dx = Amino[i].x - Amino[i+1].x;
    dy = Amino[i].y - Amino[i+1].y;
    dz = Amino[i].z - Amino[i+1].z;
    R0[i][i+1] = sqrt(dx*dx+dy*dy+dz*dz);
  }

  for(i = 1 ; i <= tot_amino ; i ++){
    Amino0[i].x = Amino[i].x;
    Amino0[i].y = Amino[i].y;
    Amino0[i].z = Amino[i].z;
    Amino0[i].a1 = Amino[i].a1;
    Amino0[i].zeta1 = Amino[i].zeta1;
    d_Amino[i].x = 0.;
    d_Amino[i].y = 0.;
    d_Amino[i].z = 0.;
  } 
  return(0);
  
}

