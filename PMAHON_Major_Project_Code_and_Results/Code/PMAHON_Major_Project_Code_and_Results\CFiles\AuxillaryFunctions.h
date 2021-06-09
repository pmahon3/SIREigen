/* Header file for AuxillaryFunctions.c which implements supplementary mathematical and statistical functions operating on the Graph structure defined in Graph.c */

double powerMethod( int** mat, int iters, int size);
int compareDouble( double a, double b);
void genBAAdjMat( int** adjmat, int size, int m0 );
void genERAdjMat( int** adjmat, int size, double p);
