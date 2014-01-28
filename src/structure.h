typedef struct Atom{

  char   name[4], chain[1];
  int    num, ind, posn; /* position and index of the whole aminoacid along the sequence */
  double x, y, z;
  double fcx, fcy, fcz;
  double a1, zeta1;

} atom;
