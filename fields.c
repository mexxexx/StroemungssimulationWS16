#include <stdio.h>
#include <stdlib.h>
#include "fields.h"
#include "matrix.h"

double *sampleFDgridOnCellCorners(double (*func)(double, double), double xlength, double ylength, int imax, int jmax) {
	double *grid = malloc(imax * jmax * sizeof(double));
	if (grid == NULL) {
		printf("Es konnte kein Speicher für das Gitter allokiert werden.\n");
		return NULL;
	}
	
	double deltaX = xlength / (imax + 1);
	double deltaY = ylength / (jmax + 1);

	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
			grid[POS2D(i, j, jmax)] = func((i+1)*deltaX, (j+1)*deltaY);
	
	return grid;
}

double *sampleFDgridOnCellCenters(double (*func)(double, double), double xlength, double ylength, int imax, int jmax) {
	double *grid = malloc((imax+2) * (jmax+2) * sizeof(double));
	if (grid == NULL) {
		printf("Es konnte kein Speicher für das Gitter allokiert werden.\n");
		return NULL;
	}
	
	double deltaX = xlength / (imax);
	double deltaY = ylength / (jmax);
	
	for (int i = 0; i < imax+2; i++) {
		for (int j = 0; j < jmax+2; j++) {
			if (i == 0 || i == imax+1 || j == 0 || j == jmax+1)
				grid[POS2D(i, j, jmax+2)] = 0;
			else
				grid[POS2D(i, j, jmax+2)] = func((i-0.5)*deltaX, (j-0.5)*deltaY);
		}
	}
	
	return grid;
}
