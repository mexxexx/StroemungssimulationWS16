#include <stdio.h>
#include <stdlib.h>
#include "poisson.h"

void freeMatrix(double** matrix, int n) {
	int i;
	for (i = n-1; i>= 0; i--)
		free(matrix[i]);
	free(matrix);
}

int getIndexInSystemMatrix(int i, int j, int jmax) {
	return i * jmax + j;
}

double** create2DpoissonMatrix(double xlength, double ylength, int imax, int jmax) {
	int matrixDim = imax * jmax;
	double** matrix = (double**)calloc(matrixDim, sizeof(double*));
	if (matrix == NULL) {
		printf ("Fehler bei Spiecherallokation für Matrix\n");
		return NULL;
	}
	
	double oneOverDeltaXSquared = ((imax + 1) * (imax + 1)) / (xlength * xlength);
	double oneOverDeltaYSquared = ((jmax + 1) * (jmax + 1)) / (ylength * ylength);
	int i, j;
	
	for (i = 0; i < matrixDim; i++) {
		matrix[i] = (double*)calloc(matrixDim, sizeof(double));
		if (matrix[i] == NULL) {
			printf ("Fehler bei Spiecherallokation für Matrixzeile %i\n", i);
			freeMatrix(matrix, i);
			return NULL;
		}
	}
	
	int rowOffset = 0;
	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax; j++) {
			double* currentRow = matrix[rowOffset++];
			currentRow[getIndexInSystemMatrix(i, j, jmax)] = 2 * (oneOverDeltaXSquared + oneOverDeltaYSquared);
			
			if (i > 0)
				currentRow[getIndexInSystemMatrix(i-1, j, jmax)] = -oneOverDeltaXSquared;
			
			if (i < imax-1)
				currentRow[getIndexInSystemMatrix(i+1, j, jmax)] = -oneOverDeltaXSquared;
			
			if (j > 0)
				currentRow[getIndexInSystemMatrix(i, j-1, jmax)] = -oneOverDeltaYSquared;
			
			if (j < jmax-1)
				currentRow[getIndexInSystemMatrix(i, j+1, jmax)] = -oneOverDeltaYSquared;
		}
	}
	
	return matrix;
}

void printMatrix(double** matrix, int n, int m) {
	int x, y;
	for (y = 0; y < n; y++) {
		for (x = 0; x < m; x++)
			printf("%f ", matrix[y][x]);
		printf("\n");
	}
}
