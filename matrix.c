#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void freeMatrix(double **matrix, int rows) {
	for (int i = rows-1; i>= 0; i--)
		free(matrix[i]);
	free(matrix);
}

void printMatrix(double *matrix, int rows, int columns) {
	for (int j = rows-1; j>=0; j--) {
		for (int i = 0; i < columns; i++) 
			printf("%f ", matrix[POS2D(i, j, columns)]);
		printf("\n");
	}
}

void printScalarField(double *matrix, int imax, int jmax, double xlength, double ylength, char *filename) {
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		return;
	}
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Scalar Field\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	fprintf(f, "DIMENSIONS %i %i 1\n", imax, jmax);
	fprintf(f, "X_COORDINATES %i double\n", imax);
	for (int i = 0; i < imax; i++) 
		fprintf(f, "%f ", i * deltaX);
	fprintf(f, "\nY_COORDINATES %i double\n", jmax);
	for (int j = 0; j < jmax; j++) 
		fprintf(f, "%f ", j * deltaY);
	fprintf(f, "\nZ_COORDINATES 1 double\n");
	fprintf(f, "0.0\n");
	fprintf(f, "POINT_DATA %i\n", imax * jmax);
	fprintf(f, "SCALARS Skalarfeld double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int j = 1; j <= jmax; j++) 
		for (int i = 1; i <= imax; i++) 
			fprintf(f, "%f\n", matrix[POS2D(i, j, imax+2)]);
	
	fclose(f);
}

void printVectorField(double *xvalues, double *yvalues, int imax, int jmax, double xlength, double ylength, char *filename) {
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		return;
	} 
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Vector Field\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	fprintf(f, "DIMENSIONS %i %i 1\n", imax, jmax);
	fprintf(f, "X_COORDINATES %i double\n", imax);
	for (int i = 0; i < imax; i++) 
		fprintf(f, "%f ", i * deltaX);
	fprintf(f, "\nY_COORDINATES %i double\n", jmax);
	for (int j = 0; j < jmax; j++) 
		fprintf(f, "%f ", j * deltaY);
	fprintf(f, "\nZ_COORDINATES 1 double\n");
	fprintf(f, "0.0\n");
	fprintf(f, "POINT_DATA %i\n", imax * jmax);
	fprintf(f, "VECTORS Vektorfeld double\n");
	for (int j = 1; j <= jmax; j++) 
		for (int i = 1; i <= imax; i++) 
			fprintf(f, "%f %f 0.0\n", xvalues[POS2D(i, j, imax+2)], yvalues[POS2D(i, j, imax+2)]);
	
	fclose(f);
}
