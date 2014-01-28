/* this is the program "write_protein.c" */

/* it writes the conformation of the chain in Rasmol format */

#include "def_param.h"


extern char Protein_name[5];

void write_protein(){

  int    i, step2;
  double min_x, min_y, min_z;
  char   filenamestr[50],folder[50];


  for(i = 1 ; i <= tot_amino ; i ++){
    nr_N[i]    = i;
    strcpy(mol_N[i],"CA");
  }
  step2=step/nav2;
  sprintf(filenamestr,"Struct_data/%d/%s.structure_%d.pdb",run,Protein_name,step2);
  fout1 = fopen(filenamestr,"w");
  
  for(i = 1 ; i <= tot_amino ; i++){
    fprintf(fout1,"ATOM  %5d  %2s %4s %s%4d  %9.3f %7.3f %7.3f %7.2f %7.2f\n", nr_N[i],mol_N[i], 
      Amino[i].name,Amino[i].chain,order_N[i],Amino[i].x,Amino[i].y,Amino[i].z,Amino[i].a1,Amino[i].zeta1); 
  }

  fclose(fout1);

}

