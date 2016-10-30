#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "poisson.h"
#include "fields.h"

/* Mit dieser Funktion wird p(0) initialisiert */
double initFunc(double x, double y) {
	return 0;
}

double testFunc(double x, double y) {
	return -(8 * M_PI * M_PI * sin(2 * M_PI * x) * sin (2 * M_PI * y));
}

double evalFunc(double x, double y) {
	return sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

int aufgabe1(double xlength, double ylength, int imax, int jmax) {
	int matSize = imax * jmax;
	
	double *A = create2DpoissonMatrix(xlength, ylength, imax, jmax);
	if (A == NULL)
		return 1;
	//printMatrix(A, matSize, matSize);
	
	double *rhsGrid = sampleFDgridOnCellCorners(testFunc, xlength, ylength, imax, jmax);
	if (rhsGrid == NULL) {
		free(A);
		return 1;
	}
		
	double *x = malloc(matSize * sizeof(double));
	double *b = malloc(matSize * sizeof(double));
	for (int j = 0; j < jmax; j++) {
		for (int i = 0; i < imax; i++) {
			int pos = POS2D(i, j, imax);
			x[pos] = 0;
			b[pos] = rhsGrid[pos];
		}
	}
	free(rhsGrid);
	
	solveSOR(A, x, b, matSize, matSize, 1, 1e-10, 10000);
	free(b);
	
	double *evalGrid = sampleFDgridOnCellCorners(evalFunc, xlength, ylength, imax, jmax);
	if (evalGrid == NULL) {
		free(A);
		free(x);
		return 1;
	}
	
	double error = 0;
	for (int j = 0; j < jmax; j++) {
		for (int i = 0; i < imax; i++) {
			int pos = POS2D(i, j, imax);
			double diff = fabs(evalGrid[pos] - x[pos]);
			error += diff * diff;
		}
	}
	
	//printMatrix(A, matSize, matSize);
	//printMatrix(x, matSize, 1);
	
	free(evalGrid);
	free(A);
	free(x);
	
	printf("Die Abweichung zwischen Lösung und Diskretisierung beträgt %f\n", sqrt(error));
	
	return 0;
}

int aufgabe2(double xlength, double ylength, int imax, int jmax) {
	double *p = sampleFDgridOnCellCenters(initFunc, xlength, ylength, imax, jmax);
	if (p == NULL)
		return 1;
	//printMatrix(p, jmax+2, imax+2);
	
	double *rhs = sampleFDgridOnCellCenters(testFunc, xlength, ylength, imax, jmax);
	if (rhs == NULL) {
		free(p);
		return 1;
	}
	
	//printMatrix(rhs, jmax+2, imax+2);
	solveSORforPoisson(p, rhs, 1.88177, 1e-8, 1e7, 0, xlength, ylength, imax, jmax);
	free(rhs);
	
	double *evalGrid = sampleFDgridOnCellCenters(evalFunc, xlength, ylength, imax, jmax);
	if (rhs == NULL)  {
		free(p);
		return 1;
	}
	
	//printMatrix(p, jmax+2, imax+2);
	//printMatrix(evalGrid, jmax+2, imax+2);
	
	double error = 0;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			double diff = fabs(evalGrid[POS2D(i, j, imax+2)] - p[POS2D(i, j, imax+2)]);
			error += diff * diff;
		}
	}
	
	free(p);
	free(evalGrid);
	
	printf("Die Abweichung zwischen Lösung und Diskretisierung beträgt %f\n", sqrt(error));
		
	return 0;
}

int main() {
	double xlength = 1;
	double ylength = 1;
	int imax;
	int jmax;
	printf("Auflösung auf der x-Achse: ");
	scanf("%i", &imax);
	printf("Auflösung auf der y-Achse: ");
	scanf("%i", &jmax);
	
	xlength = imax;
	ylength = jmax;
	//return aufgabe1(xlength, ylength, imax, jmax);
	return aufgabe2(xlength, ylength, imax, jmax);
}
