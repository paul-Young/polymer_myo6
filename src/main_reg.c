/* this is the program main.c */

/* "make" to compile "many_runs.sh" to run */
/*  edited Nov 2013 to a semi-flexible polymer diffusion */
/* input is <PDB ref file> <starting conf> <mseed value (must be negative number)>  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "def_param.h"
#include "config_reader.c"

extern void read_protein();
extern double distance();
void write_protein(), rforce(), force(), iteration();
void initialize1();

char Protein_name[15];

int main(int argc, char *argv[]){
  
	read_config("conf.txt");
	
	char   *str1, *str2, *strA, filename[100], filename_data[100],folder_name[100];
	int    i, j, pos, Ntemps, k;
	double dx, dy, dz, dr;
        double rad1;

	tempav = 0.0;

	// moved this from read_protein, must have tot_amino in def_param.h set 
	Amino   = malloc((tot_amino+1) * sizeof(atom));
	Amino0  = malloc((tot_amino+1) * sizeof(atom));
	d_Amino   = malloc((tot_amino+1) * sizeof(atom));
	order_N = (int*)malloc((tot_amino+1) * sizeof(int));

	str1 = argv[1];
	printf("1st input file is %s  \n", str1);
	read_protein(str1);
	initialize();

        // set up end tangent constraint paramaters 
        rad1=Amino[1].a1+Amino[2].a1;
        tangy=Amino[1].y;
        tangz=Amino[1].z+rad1*cos(thetac*3.141592653589793238462/180.);
        tangx=Amino[1].x+rad1*sin(thetac*3.141592653589793238462/180.);
    
        printf("the tangents are: %f, %f, %f\n", tangx, tangy, tangz);

	// read name of protein (i.e., the PDB id)   
	pos = 0;
	str1 += 3;
	while(*str1 != '.'){
		Protein_name[pos] = *str1;
		str1 ++;
		pos ++;
	}
	Protein_name[pos] = '\0';   

	str2 = argv[2];
	printf("2nd input file is %s  \n", str2);
	read_protein(str2); // reading in the 2nd structure, starting from it 

	run = atoi(argv[3]);
	mseed = run-10000; // make folders before mseed changes
	mkdir("Struct_data",0777);
	mkdir("Data",0777);
	sprintf(folder_name,"Struct_data/%d",run);
  	mkdir(folder_name,0777);
  	sprintf(filename_data,"Data/out%s_%4.2f_%d.dat",Protein_name,temp,run);
	fdata=fopen(filename_data,"w");  // open the main out_ file
        printf("... finished opening files....\n ");

	step = 0;
	printf("... step 0....\n "); 

	// turn on the random force 
	rforce(); 
	// calculate the force that acts on the protein 
	force();
	
	bound = -1; // not bound to the next site
	dr2_min=INF; // min distance to binding site
			
	while(step <= stepsim){
		iteration();
		step ++;
	}
	
	fclose(fdata);
	
	return 0;
}
