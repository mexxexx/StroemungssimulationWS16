#ifndef MATRIX_H
#define MATRIX_H

/* Rechnet eine 2D Position für den Zugriff auf ein 1D Array um */
#define POS2D(i, j, width) ((i) + ((j) * (width)))

/* Gibt eine Matrix aus */
void printMatrix(double *matrix, int rows, int columns);

/* Gibt eine zweidimensionale Matrix frei */
void freeMatrix(double **matrix, int rows);

void printScalarField(double *matrix, int imax, int jmax, double xlength, double ylength, char *filename);

void printVectorField(double *xvalues, double *yvalues, int imax, int jmax, double xlength, double ylength, char *filename);

#endif
