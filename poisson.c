#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "poisson.h"

int getIndexInSystemMatrix(int i, int j, int jmax) {
	return i * jmax + j;
}

double *create2DpoissonMatrix(double xlength, double ylength, int imax, int jmax) {
	int matrixDim = imax * jmax;
	double *matrix = (double*)calloc(matrixDim * matrixDim, sizeof(double));
	if (matrix == NULL) {
		printf ("Fehler bei Spiecherallokation für Matrix\n");
		return NULL;
	}
	
	double oneOverDeltaXSquared = ((imax + 1) * (imax + 1)) / (xlength * xlength);
	double oneOverDeltaYSquared = ((jmax + 1) * (jmax + 1)) / (ylength * ylength);
	
	/*for (i = 0; i < matrixDim; i++) {
		matrix[i] = (double*)calloc(matrixDim, sizeof(double));
		if (matrix[i] == NULL) {
			printf ("Fehler bei Spiecherallokation für Matrixzeile %i\n", i);
			freeMatrix(matrix, i);
			return NULL;
		}
	}*/
	
	int rowOffset = 0;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			//matrix[POS2D(i, j, matrixDim)] = ++rowOffset;
			double* currentRow = &matrix[rowOffset++ * matrixDim];
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

void solveSOR(double *A, double* x, double *b, int rows, int cols, double omega, double epsilon, int itermax) {
	double error = 1;
	int iter = 0;
	while (error >= epsilon && iter < itermax) {
		error = 0;
		for (int i = 0; i < cols; i++) {
			double xnew = b[i];
			for (int j = 0; j < cols; j++) {
				if (j != i)
					xnew -= A[POS2D(i, j, rows)] * x[j];
			}
			
			xnew *= omega / A[POS2D(i, i, rows)];
			xnew += (1 - omega) * x[i];
			error = fmax(error, fabs(x[i]-xnew));
			x[i] = xnew;
		}
		
		iter++;
	}
	
	if (iter == itermax) 
		printf("Abgebrochen nach %i iterationen mit einem Fehler von %f\n", iter, error);
}
