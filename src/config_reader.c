#define max 500

void read_config(char* filename){
	
	FILE *fp;
	fp = fopen(filename,"r");
	char line[max],name[max];
	double value;

	while(fgets(line,max,fp) && (line[0]!='\n')){
		sscanf(line,"%s\t%lf",name,&value);
		
		
		if (!strcmp(name,"vc")) vc=value;
		if (!strcmp(name,"thetac")) thetac=value;
		if (!strcmp(name,"pull_force")) pull_force=value;
		
	}
	
	/*
	printf("vc\t%lf\n",vc);
	printf("thetac\t%lf\n",thetac);
	printf("pull_force\t%lf\n",pull_force);
	*/

return;
}
