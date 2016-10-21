#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void freeMatrix(double **matrix, int rows) {
	for (int i = rows-1; i>= 0; i--)
		free(matrix[i]);
	free(matrix);
}

void printMatrix(double *matrix, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			printf("%f ", matrix[POS2D(i, j, rows)]);
		printf("\n");
	}
	printf("\n");
}
