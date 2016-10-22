#ifndef FIELDS_H
#define FIELD_H

/* Erstellt ein Gitter auf Eckpunkten von der gegebenen Funktion 
 * auf dem Intervall [0, xlength] x [0, ylength] mit einer Auflösung 
 * von imax Gitterpunkten in x- und jmax in y-Richtung.*/
double *sampleFDgridOnCellCorners(double (*func)(double, double), double xlength, double ylength, int imax, int jmax);

/* Erstellt ein Gitter auf Zellmittelpunkten von der gegebenen Funktion 
 * auf dem Intervall [0, xlength] x [0, ylength] mit einer Auflösung 
 * von (imax+1) Gitterpunkten in x- und (jmax+1) in y-Richtung.*/
double *sampleFDgridOnCellCenters(double (*func)(double, double), double xlength, double ylength, int imax, int jmax);

#endif
