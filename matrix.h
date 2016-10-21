#ifndef MATRIX_H
#define MATRIX_H

#define POS2D(i, j, height) (i * height + j)

void printMatrix(double *matrix, int rows, int columns);
void freeMatrix(double **matrix, int rows);

#endif
