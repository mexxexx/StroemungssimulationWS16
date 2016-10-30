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
		for (int j = 0; j < rows; j++) {
			double xnew = b[j];
			for (int i = 0; i < cols; i++) {
				if (j != i)
					xnew -= A[POS2D(i, j, cols)] * x[i];
			}
			
			xnew *= omega / A[POS2D(j, j, cols)];
			xnew += (1 - omega) * x[j];
			error = fmax(error, fabs(x[j]-xnew));
			x[j] = xnew;
		}
		
		iter++;
	}
	
	if (iter == itermax) 
		printf("Abgebrochen nach %i iterationen mit einem Fehler von %f\n", iter, error);
}

/* p ist Gitter mit (imax x jmax) inneren Zellen */
void applyHomogenousNeumannBC(double *p, int imax, int jmax) {
	const int imaxPlus2 = imax+2;
	const int imaxPlus1 = imax+1;
	const int jmaxPlus1 = jmax+1;
	const int lowerBound = (imaxPlus1 < jmaxPlus1) ? imaxPlus1 : jmaxPlus1;
	int k;
	
	for (k = 1; k < lowerBound; k++) {
		p[POS2D(k, 0, imaxPlus2)] = p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = p[POS2D(k, jmax, imaxPlus2)];
		p[POS2D(0, k, imaxPlus2)] = p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = p[POS2D(imax, k, imaxPlus2)];
	}
	
	for (k = lowerBound; k < imaxPlus1; k++) {
		p[POS2D(k, 0, imaxPlus2)] = p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = p[POS2D(k, jmax, imaxPlus2)];
	}
	
	for (k = lowerBound; k < jmaxPlus1; k++) {
		p[POS2D(0, k, imaxPlus2)] = p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = p[POS2D(imax, k, imaxPlus2)];
	}
}

void applyHomogenousDirichletBC(double *p, int imax, int jmax) {
	const int imaxPlus2 = imax+2;
	const int imaxPlus1 = imax+1;
	const int jmaxPlus1 = jmax+1;
	const int lowerBound = (imaxPlus1 < jmaxPlus1) ? imaxPlus1 : jmaxPlus1;
	int k;
	
	for (k = 1; k < lowerBound; k++) {
		p[POS2D(k, 0, imaxPlus2)] = -p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = -p[POS2D(k, jmax, imaxPlus2)];
		p[POS2D(0, k, imaxPlus2)] = -p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = -p[POS2D(imax, k, imaxPlus2)];
	}
	
	for (k = lowerBound; k < imaxPlus1; k++) {
		p[POS2D(k, 0, imaxPlus2)] = -p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = -p[POS2D(k, jmax, imaxPlus2)];
	}
	
	for (k = lowerBound; k < jmaxPlus1; k++) {
		p[POS2D(0, k, imaxPlus2)] = -p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = -p[POS2D(imax, k, imaxPlus2)];
	}
}

void solveSORforPoisson(double *p, double *rhs, double omega, double epsilon, int itermax, 
						int useNeumannBC, double xlength, double ylength, int imax, int jmax) {
	const double oneMinusOmega = 1 - omega;
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	const double oneOverDeltaXSquared = 1 / (deltaX * deltaX);
	const double oneOverDeltaYSquared = 1 / (deltaY * deltaY);
	const double omegaRelaxation = omega / (2 * (oneOverDeltaXSquared + oneOverDeltaYSquared));
	const int imaxPlus2 = imax+2;
	
	int ijlocation = 0;
	
	double error = 1;
	double sum = 0;
	double twoPij = 0;
	int iter = 0;
	int i,j;
	while (error >= epsilon && iter < itermax) {
		error = 0;
		if (useNeumannBC)
			applyHomogenousNeumannBC(p, imax, jmax);
		else
			applyHomogenousDirichletBC(p, imax, jmax);
			
		for (j = 1; j <= jmax; j++) {
			for (i = 1; i <= imax; i++) {				
				ijlocation = POS2D(i, j, imaxPlus2);
				p[ijlocation] = oneMinusOmega * p[ijlocation] + omegaRelaxation * (
						oneOverDeltaXSquared * (p[ijlocation-1] + p[ijlocation+1]) + 
						oneOverDeltaYSquared * (p[ijlocation-imaxPlus2] + p[ijlocation+imaxPlus2]) - 
						rhs[ijlocation]);
			}
		}
		
		for (i = 1; i <= imax; i++) {
			for (j = 1; j <= jmax; j++) {
				ijlocation = POS2D(i, j, imaxPlus2);
				twoPij = 2 * p[ijlocation];
				sum = oneOverDeltaXSquared * (p[ijlocation+1] + p[ijlocation-1] - twoPij) +
					oneOverDeltaYSquared * (p[ijlocation+imaxPlus2] + p[ijlocation-imaxPlus2] - twoPij) - 
					rhs[ijlocation];
				error += sum * sum;
			}
		}
		error /= imax * jmax;
		error = sqrt(error);
		
		iter++;
	}
	
	if (iter == itermax) 
		printf("Abgebrochen nach %i iterationen mit einem Fehler von %f\n", iter, error);
}
