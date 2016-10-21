#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "poisson.h"
#include "fields.h"

double testFunc(double x, double y) {
	return 8 * M_PI * M_PI * sin(2 * M_PI * x) * sin (2 * M_PI * y);
}

double evalFunc(double x, double y) {
	return sin(2 * M_PI * x) * sin(2 * M_PI * y);
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
	
	int matSize = imax * jmax;
	
	double *A = create2DpoissonMatrix(xlength, ylength, imax, jmax);
	if (A == NULL)
		return 1;
		
	double *testGrid = sampleFDgridOnCellCorners(testFunc, xlength, ylength, imax, jmax);
	if (testGrid == NULL)
		return 1;
	//printMatrix(testGrid, imax, jmax);
		
	double *x = malloc(matSize * sizeof(double));
	double *b = malloc(matSize * sizeof(double));
	int curr = 0;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			x[curr] = 0;
			b[curr] = testGrid[POS2D(i, j, jmax)];
			curr++;
		}
	}
	free(testGrid);
	
	solveSOR(A, x, b, matSize, matSize, 1, 1e-10, 10000);
	
	double *evalGrid = sampleFDgridOnCellCorners(evalFunc, xlength, ylength, imax, jmax);
	curr = 0;
	double error = 0;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			double diff = fabs(evalGrid[POS2D(i, j, jmax)] - x[curr]);
			error += diff * diff;
			curr++;
		}
	}
	free(evalGrid);
	
	printf("Die Abweichung zwischen Lösung und Diskretisierung beträgt %f\n", sqrt(error));
	
	//printMatrix(A, matSize, matSize);
	//printMatrix(x, 1, matSize);
	
	free(A);
	free(x);
	free(b);
	
	
	return 0;
}
