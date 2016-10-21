#ifndef POISSON_H
#define POISSON_H

void freeMatrix(double** matrix, int n);
int getIndexInSystemMatrix(int i, int j, int jmax);
double** create2DpoissonMatrix(double ilength, double jlength, int imax, int jmax);
void printMatrix(double** matrix, int n, int m);

#endif
