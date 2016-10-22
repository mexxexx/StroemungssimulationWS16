#ifndef MATRIX_H
#define MATRIX_H

/* Rechnet eine 2D Position für den Zugriff auf ein 1D Array um */
#define POS2D(i, j, width) (((width) * (i)) + j)

/* Gibt eine Matrix aus */
void printMatrix(double *matrix, int rows, int columns);

/* Gibt eine zweidimensionale Matrix frei */
void freeMatrix(double **matrix, int rows);

#endif
