#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>

#include "def_param.h"
#define max 500

extern int     string_to_integer();

int read_protein(char *filename_PRO){

  int    i,j, k, pos;
  char   buffer[max], *s;
  FILE   *fp;
  double rg, rgsq;


  fp = fopen(filename_PRO,"r");  
  // Check # of atoms
  i = 0;
  while(fgets(buffer,max,fp) && (buffer[0]!='\n')){
    s = buffer;
    if(!strncmp(s,"ATOM",4)){
      s += 13;
      if(!strncmp(s,"CA",2)) i ++; 
    }
    if(!strncmp(s,"TER",3)) break; 
  }

  if (i != tot_amino ){
    printf("wrong number of amino acids, tot_amino = %d, i = %d, exiting \n", tot_amino, i);
    exit(0);
  }
 
  rewind(fp);
  
  // read file
  for(i = 1 ; i <= tot_amino; i++){
    Amino[i].x = 0.0;
    Amino[i].y = 0.0;
    Amino[i].z = 0.0;
  }
 
  i = 0;
  while(fgets(buffer,max,fp)&& (buffer[0]!='\n')){
    s = buffer; 
    if(!strncmp(s,"ATOM",4)){
      s += 13;
      if(!strncmp(s,"CA",2)){ 
        s += 4;
        i ++;    // start counting from 1 
        sscanf(s,"%3s",Amino[i].name);
        //printf("res name: %3s \n", Amino[i].name);
        s += 4;
        sscanf(s,"%s",Amino[i].chain);
        //printf("%s",s);
        //printf("chain is %s \n", Amino[i].chain);
        sscanf(s,"%*s %d",&order_N[i]);
        s += 5;
        sscanf(s,"%lf %lf %lf",&Amino[i].x,&Amino[i].y,&Amino[i].z);
        //printf("%f %f %f \n", Amino[i].x,Amino[i].y,Amino[i].z);
        s += 5;
        sscanf(s,"%lf %lf %lf",&Amino[i].x,&Amino[i].y,&Amino[i].z);
        //printf("%f %f %f \n", Amino[i].x,Amino[i].y,Amino[i].z);
        s += 24;
        //printf("%s",s);
        sscanf(s,"%lf %lf",&Amino[i].a1,&Amino[i].zeta1);
        printf("%f %f \n", Amino[i].a1,Amino[i].zeta1);
        
      }  
    }
    if(!strncmp(s,"TER",3)) break;
  }
  fclose(fp);
 
  /* get the square of the radius of gyration of the protein */
  rgsq = 0.0; 
  for(i = 1 ; i <= (tot_amino-1) ; i ++){
    for(j = (i+1) ; j <= tot_amino ; j ++){

      rgsq = rgsq + (Amino[i].x-Amino[j].x)*(Amino[i].x-Amino[j].x) + 
	(Amino[i].y-Amino[j].y)*(Amino[i].y-Amino[j].y) + (Amino[i].z-Amino[j].z)*(Amino[i].z-Amino[j].z);
    }
  }

  rgsq = rgsq/((double)(tot_amino*tot_amino));
  rg = sqrt(rgsq);
  printf("rg = %f\n", rg);

  return(0);

}
