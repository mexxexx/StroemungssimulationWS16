#ifndef POISSON_H
#define POISSON_H

/* Gibt die Position eines Wertes (i,j) in der aktuellen Zeile der Systemmatrix zurück */
int getIndexInSystemMatrix(int i, int j, int jmax);

/* Erstellt eine Poissonmatrix auf dem Gebiet [0, xlength] x [0, ylength] 
 * mit einer Auflösung von imax Gitterpunkten in x- und jmax in y-Richtung.*/
double *create2DpoissonMatrix(double xlength, double ylength, int imax, int jmax);

/* Löst das LGS Ax=b mit dem SOR Verfahren */
void solveSOR(double *A, double *x, double *b, int rows, int cols, double omega, double epsilon, int itermax);

/* Wendet die homogene Neumann-Randbedingung Laplace(p|Gamma)=0 
 * auf eine Poissonmatrix p mit (imax x jmax) inneren Gitterpunkten mit Ghost-Schicht an. */
void applyHomogenousNeumannBC(double *p, int imax, int jmax);

/* Wendet die homogene Dirichlet-Randbedingung p|Gamma=0 
 * auf eine Poissonmatrix p mit (imax x jmax) inneren Gitterpunkten mit Ghost-Schicht an. */
void applyHomogenousDirichletBC(double *p, int imax, int jmax);

/* Löst das LGS Laplace(p)=rhs mit einem effizienten SOR Verfahren, 
 * wobei p eine (imax+2) x (jmax+2) Matrix ist die imax bzw. jmax 
 * innere Gitterpunkte und eine Ghost-Schicht besitzt. */
void solveSORforPoisson(double *p, double *rhs, double omega, double epsilon, int itermax, 
						int useNeumannBC, double xlength, double ylength, int imax, int jmax);

#endif
