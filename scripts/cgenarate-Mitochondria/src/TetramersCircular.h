#include "settings.h"
extern int nthreads;

#include <stdio.h>

#ifndef TETRAMERS_CIRCULAR_H
#define TETRAMERS_CIRCULAR_H

void ToTetramersCircular(char *RNAME, int n, int *TNums);
void MassFromResidueCircular(char *RNAME, int n, double *mass);

void LinkerForcesTetramersCircular(double (*F)[3], double (*r)[3], int n, int nn, int *TNums,
								   double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
								   double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
								   double KLong[4][Order]);
void EnergyTetramersCircular(double (*r)[3], int n, int nn, int *TNums,
							 double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
							 double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
							 double KLong[4][Order], double *const E);

void PrintPlotsCircular(double (*r)[3], int n, int nn, FILE *plotfile);

#endif // TETRAMERS_CIRCULAR_H
