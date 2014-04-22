#ifndef PROTEIN_H
#define PROTEIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ostream>

#define BUF_MAX 100

class Amino
{
public:
    char    name[4],chain;
    int     idx, posn;
    double  x,y,z,fcx,fcy,fcz,a1,zeta1;

    Amino();
};

class Protein
{
public:
    std::string name;
    int size;
    Amino *chain;

    Protein();
    Protein(std::string);
    void read(std::string);
    double distance(int,int);
};

std::ostream& operator<<(std::ostream& os, const Protein& obj);
std::ostream& operator<<(std::ostream& os, const Amino& obj);

#endif // PROTEIN_H


