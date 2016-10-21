#include <stdio.h>
#include <stdlib.h>
#include "fields.h"
#include "matrix.h"

double *sampleFDgridOnCellCorners(double (*func)(double, double), double xlength, double ylength, int imax, int jmax) {
	double *grid = malloc(imax * jmax * sizeof(double*));
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
