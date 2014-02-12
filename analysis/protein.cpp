#include "protein.h"

/* -------------- Amino Constructor -------------- */
Amino::Amino()
{
    idx=-1;
    posn=-1;
    x=0;y=0;z=0;fcx=0;fcy=0;fcz=0;
}

/* -------------- Protein Constructor -------------- */
Protein::Protein()
{
}

Protein::Protein(std::string filename)
{
    sscanf(name,"myoVI");
    this->read(filename);
}

/* -------------- Protein Methods -------------- */

void Protein::read(std::string filename){
    FILE   *fp;
    char   buffer[BUF_MAX], *s;
    fp = fopen(filename.c_str(),"r");

    size=0;
    while(fgets(buffer,BUF_MAX,fp) && (buffer[0]!='\n')){
        s = buffer;
        if(!strncmp(s,"ATOM",4)){
          s += 13;
          if(!strncmp(s,"CA",2)) size++;
        }
        if(!strncmp(s,"TER",3)) break;
    }

    chain = (Amino*)malloc((size+1)*sizeof(Amino)); // start counting from 1

    rewind(fp);
    int i=0;
    while(fgets(buffer,BUF_MAX,fp)&& (buffer[0]!='\n')){
        s = buffer;
        if(!strncmp(s,"ATOM",4)){
          s += 13;
          if(!strncmp(s,"CA",2)){
            s += 4;
            i ++;    // start counting from 1
            sscanf(s,"%3s",chain[i].name);
            s += 4;
            sscanf(s,"%s",&chain[i].chain);
            s += 10;
            sscanf(s,"%lf %lf %lf",&chain[i].x,&chain[i].y,&chain[i].z);
            s += 24;
            sscanf(s,"%lf %lf",&chain[i].a1,&chain[i].zeta1);
          }
        }
        if(!strncmp(s,"TER",3)) break;
      }
      fclose(fp);

}

double Protein::distance(int i, int j){
    double d2=0;
    d2=(chain[i].x-chain[j].x)*(chain[i].x-chain[j].x)
            +(chain[i].y-chain[j].y)*(chain[i].y-chain[j].y)
            +(chain[i].z-chain[j].z)*(chain[i].z-chain[j].z);
    return sqrt(d2);
}

std::ostream& operator<<(std::ostream& os, const Protein& obj)
{
    for (int i=1;i<=obj.size;i++){
        os << obj.chain[i].z;
    }
    return os;
}
