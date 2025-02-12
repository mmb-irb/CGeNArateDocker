#include "settings.h"
#include <stdio.h>
extern int nthreads;

#ifndef TETRAMERS_LINEAR_H
#define TETRAMERS_LINEAR_H

void ToTetramersLinear(char *RNAME, int n, int *TNums);
void MassFromResidueLinear(char *RNAME, int n, double *mass);

void LinkerForcesTetramersLinear(double (*F)[3], double (*r)[3], int n, int nn, int *TNums,
								 double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
								 double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
								 double KLong[4][Order]);
void EnergyTetramersLinear(double (*r)[3], int n, int nn, int *TNums,
						   double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
						   double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
						   double KLong[4][Order], double *const E);

void PrintPlotsLinear(double (*r)[3], int n, int nn, FILE *plotfile);

#endif // TETRAMERS_LINEAR_H
