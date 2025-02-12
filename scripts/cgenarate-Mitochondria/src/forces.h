#include "settings.h"
extern int nthreads;

#ifndef FORCES_H
#define FORCES_H

/*Force Functions*/
void BondForceij(double (*F)[3], double (*r)[3], double d, double *K, int i, int j);
void BondForce(double F[3], double Fs[3], double r[3], double s[3], double d, double K[Order]);
void AngleForceijk(double (*F)[3], double (*r)[3], double angle, double *K, int i, int j, int k);
void AngleForce(double F0[3], double F1[3], double F2[3], double r0[3], double r1[3], double r2[3], double angle,
				double K[Order]);
void DH(double (*F)[3], double (*r)[3], double q, double Kappa, double Eps, int n);
void LJ(double (*F)[3], double (*r)[3], double K, double sigma, int n);

void Brownian(double sigma, double *mass, double (*F)[3], int n);

/*Energy Functions*/
double BondEnergyij(double (*r)[3], double d, double *K, int i, int j);
double BondEnergy(double r[3], double s[3], double d, double K[Order]);
double AngleEnergyijk(double (*r)[3], double angle, double *K, int i, int j, int k);
double AngleEnergy(double r0[3], double r1[3], double r2[3], double angle, double K[Order]);
void DHEnergy(double (*r)[3], double q, double IonicStrength, double Eps, int n, double *E);
double LJEnergy(double (*r)[3], double EpsLJ, double sigma, int n);

double kineticEnergy(double (*v)[3], int n, double *m);

#endif // FORCES_H
