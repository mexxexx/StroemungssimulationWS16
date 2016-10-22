#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "poisson.h"
#include "fields.h"

/* Mit dieser Funktion wird p(0) initialisiert */
double initFunc(double x, double y) {
	return 1;
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
		
	double *rhsGrid = sampleFDgridOnCellCorners(testFunc, xlength, ylength, imax, jmax);
	if (rhsGrid == NULL) {
		free(A);
		return 1;
	}
	//printMatrix(testGrid, imax, jmax);
		
	double *x = malloc(matSize * sizeof(double));
	double *b = malloc(matSize * sizeof(double));
	int curr = 0;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			x[curr] = 0;
			b[curr] = rhsGrid[POS2D(i, j, jmax)];
			curr++;
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
	curr = 0;
	double error = 0;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			double diff = fabs(evalGrid[POS2D(i, j, jmax)] - x[curr]);
			error += diff * diff;
			curr++;
		}
	}
	
	//printMatrix(A, matSize, matSize);
	//printMatrix(x, 1, matSize);
	
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
	
	double *rhs = sampleFDgridOnCellCenters(testFunc, xlength, ylength, imax, jmax);
	if (rhs == NULL) {
		free(p);
		return 1;
	}
	
	solveSORforPoisson(p, rhs, 1, 1e-10, 1e7, 0, xlength, ylength, imax, jmax);
	free(rhs);
	
	double *evalGrid = sampleFDgridOnCellCenters(evalFunc, xlength, ylength, imax, jmax);
	if (rhs == NULL)  {
		free(p);
		return 1;
	}
	
	//printMatrix(p, imax +2, jmax+2);
	//printMatrix(evalGrid, imax +2, jmax+2);
	
	double error = 0;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			double diff = fabs(evalGrid[POS2D(i, j, jmax+2)] - p[POS2D(i, j, jmax+2)]);
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
	
	//return aufgabe1(xlength, ylength, imax, jmax);
	
	return aufgabe2(xlength, ylength, imax, jmax);
}
