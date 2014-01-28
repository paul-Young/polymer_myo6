/* this is the program "force.c": it calculates the force that acts on the chain at each moment in time */

#include "def_param.h"

void bond(),repel(),bend(),bend3(),end_tangent(),pull();

void force(){
	int i,j;
	epot = 0.0;
		
	repel();
	bond();
	epot_bend=0.0;
	bend3(1,tot_pol1);
	bend3(tot_pol1+1,tot_amino);
	end_tangent();
	pull(); // apply "pull_force" in z direction for bead "tot_pol1" 

}

void pull(){
	
	Amino[tot_pol1].fcz+=pull_force;
	
}

/* get the repulsive contribution to avoid overlaps */
void repel(){

  double dx, dy, dz, rsqi, rsqi6, p6=0.0, pot, invrsqi ;
  double dphi, fx, fy, fz;
  int    i, j;

  epot_rep = 0.0;

	for(i = 1 ; i <= (tot_amino-2) ; i ++){
		for(j = (i+2) ; j <= tot_amino ; j ++){
		  dx = Amino[j].x - Amino[i].x;
		  dy = Amino[j].y - Amino[i].y;
		  dz = Amino[j].z - Amino[i].z;                  
     
		  rsqi   = (R0[i][j] * R0[i][j])/(dx*dx+dy*dy+dz*dz);
        //printf("in repel1: AA %d-%d rsqi (should be close to 1) is %f\n", i,j,rsqi);
        		
		  invrsqi  = 1.0/rsqi;
		  if( invrsqi > R_sigma) continue;				
	     rsqi6  = rsqi*rsqi*rsqi;
	     p6  = el * rsqi6;  
		  pot = p6;
		  epot_rep = epot_rep + pot;
		  dphi = rsqi * 6.0 * p6/(R0[i][j]*R0[i][j]);
	          
        fx = dx * dphi;
        fy = dy * dphi;
        fz = dz * dphi;

        Amino[i].fcx = Amino[i].fcx - fx;
        Amino[i].fcy = Amino[i].fcy - fy;
        Amino[i].fcz = Amino[i].fcz - fz;
        Amino[j].fcx = Amino[j].fcx + fx;
        Amino[j].fcy = Amino[j].fcy + fy;
        Amino[j].fcz = Amino[j].fcz + fz;    

        epot_rep = epot_rep + pot;
      
        //printf("in repel2: AA %d-%d rsqi (should be close to 1) is %f\n", i,j,rsqi);
        
      }
    }
  
}

/* get the bonded potential */
void bond(){

  int    i;
  double dx,dy,dz,pot,fx,fy,fz,dr,dR;

  epot_bond = 0.0;
 
  for(i = 1 ; i <= (tot_amino-1) ; i ++){

      dx = Amino[i+1].x - Amino[i].x;
      dy = Amino[i+1].y - Amino[i].y;
      dz = Amino[i+1].z - Amino[i].z;
      dr = sqrt(dx*dx+dy*dy+dz*dz);
      dR = (dr - R0[i][i+1]);
 
      pot = (kspring_cov/2.0) * dR * dR;
      
      if(dR > R_limit){
        printf("IN BOND\n");
        printf("%lld %d %f %f %f %f %f %f %f %f %f\n", step, i, dR, R_limit, pot,
	       d_Amino[i].x, d_Amino[i].y, d_Amino[i].z, d_Amino[i+1].x, d_Amino[i+1].y, d_Amino[i+1].z);
        printf("the distance is %f \n", R0[i][i+1]);
        exit(0);
      }
                 
      epot_bond = epot_bond + pot;

      /* the corresponding force */
      fx = kspring_cov * (1. - (R0[i][i+1])/dr) * dx;
      fy = kspring_cov * (1. - (R0[i][i+1])/dr) * dy;
      fz = kspring_cov * (1. - (R0[i][i+1])/dr) * dz;
      
      Amino[i].fcx = Amino[i].fcx + fx;
      Amino[i].fcy = Amino[i].fcy + fy;
      Amino[i].fcz = Amino[i].fcz + fz;
 
      Amino[i+1].fcx = Amino[i+1].fcx - fx;
      Amino[i+1].fcy = Amino[i+1].fcy - fy;
      Amino[i+1].fcz = Amino[i+1].fcz - fz;            
    
    
      //printf("in bond: AA %d-%d bond should be %f and is %f\n", i,i+1,R0[i][i+1],dr);
  }
}

void bend(){
	/* stiff polymer bending contraint, must have at least 5 beads to make sense */
	if (tot_amino<5){
		printf("Only have %d beads, bending constraint can't be applied!(in force_1state.c - void bend())\n",tot_amino);
		exit(1);
	}
	
	epot_bend=0.0;
	int i;
	i=1;
	Amino[i].fcx+=k_bend*(Amino[i+1].x-Amino[i+2].x);
	Amino[i].fcy+=k_bend*(Amino[i+1].y-Amino[i+2].y);
	Amino[i].fcz+=k_bend*(Amino[i+1].z-Amino[i+2].z);
	i=2;
	Amino[i].fcx+=k_bend*(Amino[i-1].x-2*Amino[i].x+2*Amino[i+1].x-Amino[i+2].x);
	Amino[i].fcy+=k_bend*(Amino[i-1].y-2*Amino[i].y+2*Amino[i+1].y-Amino[i+2].y);
	Amino[i].fcz+=k_bend*(Amino[i-1].z-2*Amino[i].z+2*Amino[i+1].z-Amino[i+2].z);
	
	for(i=3;i<=(tot_amino-2);i++){
	  Amino[i].fcx+=k_bend*(-Amino[i-2].x+2*Amino[i-1].x-2*Amino[i].x+2*Amino[i+1].x-Amino[i+2].x);
	  Amino[i].fcy+=k_bend*(-Amino[i-2].y+2*Amino[i-1].y-2*Amino[i].y+2*Amino[i+1].y-Amino[i+2].y);
	  Amino[i].fcz+=k_bend*(-Amino[i-2].z+2*Amino[i-1].z-2*Amino[i].z+2*Amino[i+1].z-Amino[i+2].z);
	}
	i=tot_amino-1;
	Amino[i].fcx+=k_bend*(-Amino[i-2].x+2*Amino[i-1].x-2*Amino[i].x+Amino[i+1].x);
	Amino[i].fcy+=k_bend*(-Amino[i-2].y+2*Amino[i-1].y-2*Amino[i].y+Amino[i+1].y);
	Amino[i].fcz+=k_bend*(-Amino[i-2].z+2*Amino[i-1].z-2*Amino[i].z+Amino[i+1].z);
	i=tot_amino;
	Amino[i].fcx+=k_bend*(-Amino[i-2].x+Amino[i-1].x);
	Amino[i].fcy+=k_bend*(-Amino[i-2].y+Amino[i-1].y);
	Amino[i].fcz+=k_bend*(-Amino[i-2].z+Amino[i-1].z);


}

void bend3(int istart, int iend){
	/* stiff polymer bending contraint, must have at least 5 beads to make sense */
	if ((iend-istart)< 3){
		printf("Only have %d-%d beads, bending constraint can't be applied!(in force_1state.c - void bend())\n",iend,istart);
		exit(1);
	}
	
	int i;

	i=istart;
	Amino[i].fcx+=k_bend*(Amino[i+1].x-Amino[i+2].x);
	Amino[i].fcy+=k_bend*(Amino[i+1].y-Amino[i+2].y);
	Amino[i].fcz+=k_bend*(Amino[i+1].z-Amino[i+2].z);
	i=istart+1;
	Amino[i].fcx+=k_bend*(Amino[i-1].x-2*Amino[i].x+2*Amino[i+1].x-Amino[i+2].x);
	Amino[i].fcy+=k_bend*(Amino[i-1].y-2*Amino[i].y+2*Amino[i+1].y-Amino[i+2].y);
	Amino[i].fcz+=k_bend*(Amino[i-1].z-2*Amino[i].z+2*Amino[i+1].z-Amino[i+2].z);
	for(i=(istart+2);i<=(iend-2);i++){
	  Amino[i].fcx+=k_bend*(-Amino[i-2].x+2*Amino[i-1].x-2*Amino[i].x+2*Amino[i+1].x-Amino[i+2].x);
	  Amino[i].fcy+=k_bend*(-Amino[i-2].y+2*Amino[i-1].y-2*Amino[i].y+2*Amino[i+1].y-Amino[i+2].y);
	  Amino[i].fcz+=k_bend*(-Amino[i-2].z+2*Amino[i-1].z-2*Amino[i].z+2*Amino[i+1].z-Amino[i+2].z);
	}
	i=iend-1;
	Amino[i].fcx+=k_bend*(-Amino[i-2].x+2*Amino[i-1].x-2*Amino[i].x+Amino[i+1].x);
	Amino[i].fcy+=k_bend*(-Amino[i-2].y+2*Amino[i-1].y-2*Amino[i].y+Amino[i+1].y);
	Amino[i].fcz+=k_bend*(-Amino[i-2].z+2*Amino[i-1].z-2*Amino[i].z+Amino[i+1].z);
	i=iend;
	Amino[i].fcx+=k_bend*(-Amino[i-2].x+Amino[i-1].x);
	Amino[i].fcy+=k_bend*(-Amino[i-2].y+Amino[i-1].y);
	Amino[i].fcz+=k_bend*(-Amino[i-2].z+Amino[i-1].z);

}

void end_tangent(){

	epot_tangent=0.;
	int i=2;
//	epot_tangent=v*((Amino01[i].x-Amino01[i-1].x)*(Amino[i].x-Amino[i-1].x)+(Amino01[i].y-Amino01[i-1].y)*(Amino[i].y-Amino[i-1].y)+(Amino01[i].z-Amino01[i-1].z)*(Amino[i].z-Amino[i-1].z)); // E=-v(x0x+y0y+z0z)
// tangx, etc are the "ideal" locations for the 2nd bead given a preferential angle
	Amino[i].fcx+=vc*(tangx-Amino[i].x);
	Amino[i].fcy+=vc*(tangy-Amino[i].y);
	Amino[i].fcz+=vc*(tangz-Amino[i].z);
}

void check(){
	double dr2;
	dr2=(xb-Amino[tot_amino].x)*(xb-Amino[tot_amino].x)+(yb-Amino[tot_amino].y)*(yb-Amino[tot_amino].y)+(zb-Amino[tot_amino].z)*(zb-Amino[tot_amino].z);
	if (dr2<DR2_limit){
		bound=step;
                update();
                write_protein();
                printf("STEPPED!\n");
		exit(0);
	}
}
