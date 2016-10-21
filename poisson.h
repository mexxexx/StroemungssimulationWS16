#ifndef POISSON_H
#define POISSON_H

int getIndexInSystemMatrix(int i, int j, int jmax);
double *create2DpoissonMatrix(double ilength, double jlength, int imax, int jmax);
void solveSOR(double *A, double* x, double *b, int rows, int cols, double omega, double epsilon, int itermax);

#endif
