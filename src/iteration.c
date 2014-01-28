
/* this is the program "iteration.c" */

#include "def_param.h"


extern void rforce(),force(),update(),check();

void iteration(){  

  int i;

  for(i = 2 ; i <= tot_amino ; i ++){  
    d_Amino[i].x = (h * Amino[i].fcx)/Amino0[i].zeta1;
    d_Amino[i].y = (h * Amino[i].fcy)/Amino0[i].zeta1;
    d_Amino[i].z = (h * Amino[i].fcz)/Amino0[i].zeta1;

    /* change of position for each bead at each time step due to the action of force */
    Amino[i].x = Amino[i].x + d_Amino[i].x;
    Amino[i].y = Amino[i].y + d_Amino[i].y;
    Amino[i].z = Amino[i].z + d_Amino[i].z;    
  }

  for(i = 1 ; i <= tot_amino ; i ++){  
    /* scale back to the first AA location:   */
  //  Amino[i].x = Amino[i].x - Amino[1].x;
  //  Amino[i].y = Amino[i].y - Amino[1].y;
  //  Amino[i].z = Amino[i].z - Amino[1].z;    
  }
  
  //printf("in iteration, step is %lld  \n", step);
  rforce();
  force();
  update();
  if (bound==-1) check();
  //printf("done with forces and update, step is  %lld  \n", step);
  

}
