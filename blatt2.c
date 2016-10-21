#include "poisson.h"

int main() {
	double xlength = 0.5;
	double ylength = 1.5;
	int imax = 3;
	int jmax = 2;
	
	double** matrix = create2DpoissonMatrix(xlength, ylength, imax, jmax);
	printMatrix(matrix, imax * jmax, imax * jmax);
	freeMatrix(matrix, imax * jmax);
	return 0;
}
