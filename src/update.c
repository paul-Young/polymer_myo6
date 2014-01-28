/* this is the program "update.c" */

#include "def_param.h"

extern void  tred2(), tqli();

extern void write_protein();
double distance();
double R;


void update(){

  int    i, j, count, count2;
  double dpossum, dx, dy, dz, dr, R1, R2, dZ;
  double E_t=0., E_p=0., E_k=0.;


  dpossum = 0.0;
  for(i = 1 ; i <= tot_amino ; i ++){
    dpossum = dpossum + d_Amino[i].x*d_Amino[i].x + d_Amino[i].y*d_Amino[i].y + d_Amino[i].z*d_Amino[i].z;
  }
  dpossum = dpossum*zeta/(2*h);
  tempav  = tempav + dpossum/(3.0*((double)(tot_amino)));   // keeping track of the temperature
 
  /* every "nav" number of steps update the number of native non-bonded contacts left and output
     the energies */
  if(step%nav1 == 0){
     
     /* get the value of the current end-to-end distance */
     R = distance(Amino[tot_amino].x,Amino[tot_amino].y,Amino[tot_amino].z,Amino[1].x,Amino[1].y,Amino[1].z);
     R1 = distance(Amino[1].x,Amino[1].y,Amino[1].z,Amino[tot_pol1].x,Amino[tot_pol1].y,Amino[tot_pol1].z);
     R2 = distance(Amino[tot_amino].x,Amino[tot_amino].y,Amino[tot_amino].z,Amino[tot_pol1].x,Amino[tot_pol1].y,Amino[tot_pol1].z);
     dr2_min = dr2_min<distance(Amino[tot_amino].x,Amino[tot_amino].y,Amino[tot_amino].z,xb,yb,zb)?dr2_min:distance(Amino[tot_amino].x,Amino[tot_amino].y,Amino[tot_amino].z,xb,yb,zb);
     dZ = Amino[tot_amino].z-Amino[1].z;

     tempav = tempav/((double)nav);
     fflush(stdout);

     E_p=epot_rep+epot_bond+epot_bend;
     E_k=3./2.*tot_amino*tempav;
     E_t=E_p+E_k;
    
     fprintf(fdata, " %lld %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %6.3f %5.3f %5.3f \n", 
         step,tempav,R,R1,R2,epot_rep,epot_bond,epot_bend,dr2_min,dZ);

     if(step%nav2 == 0){
       printf("calling write step is %lld(%3.2f%% complete), nav2 is %d \n", step, 100*((double)step)/stepsim,nav2); 
       write_protein();
      
     }  

    tempav = 0.0;
  
  }/* end step number is a multiple of "nav" */

}

double distance(double x1, double y1, double z1, double x2,double y2,double z2){
  
  double r, dx, dy, dz;
  
  dx = x1 - x2;
  dy = y1 - y2;
  dz = z1 - z2;
  r = sqrt(dx*dx+dy*dy+dz*dz);

  return(r);
}

