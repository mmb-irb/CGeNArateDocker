#include "forces.h"
#include "linear_algebra.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>

void BondForceij(double (*F)[3], double (*r)[3], double d, double *K, int i, int j)
{
	/*Applies BondForce to the corresponding Points*/
	BondForce(F[i], F[j], r[i], r[j], d, K);
}

void BondForce(double F[3], double Fs[3], double r[3], double s[3], double d, double K[Order])
{
	/*Computes the Forces associated with a Bonded interaction between 2 particles,
	with equilibrium distance d and polynomial Constants K*/
	double aux, aux2;
	int l, i;

	aux = dist(r, s);
	if (aux < tol)
		return;

	for (i = Order, aux2 = 0.; i > 0; i--)
		aux2 = (aux2 + K[i - 1] * (i + 1.)) * (aux - d);
	aux2 /= aux;
	for (l = 0; l < 3; l++)
	{
		aux = aux2 * (s[l] - r[l]);
		F[l] += aux;
		Fs[l] -= aux;
	}
}

void AngleForceijk(double (*F)[3], double (*r)[3], double angle, double *K, int i, int j, int k)
{
	/*Applies AngleForce to the corresponding Points*/
	AngleForce(F[i], F[j], F[k], r[i], r[j], r[k], angle, K);
}

void AngleForce(double F0[3], double F1[3], double F2[3], double r0[3], double r1[3], double r2[3], double angle,
				double K[Order])
{
	/*Computes the Forces associated with an Angle interaction between 3 particles,
	with equilibrium angle 'angle' and polynomial Constants K*/
	double v[3], vp[3];
	double normscal, normv, normvp, frac;
	int l, i;

	for (l = 0; l < 3; l++)
	{
		v[l] = r0[l] - r1[l];
		vp[l] = r2[l] - r1[l];
	}

	normv = norm(v);
	normvp = norm(vp);

	for (l = 0; l < 3; l++)
	{
		v[l] /= normv;
		vp[l] /= normvp;
	}

	normscal = scalar(v, vp);

	for (i = Order, frac = 0.; i > 0; i--)
		frac = (frac + K[i - 1] * (i + 1.)) * (acos(normscal) - angle);

	frac /= sqrt(1 - normscal * normscal);

	for (l = 0; l < 3; l++)
	{
		F0[l] += frac * (vp[l] - normscal * v[l]) / normv;
		F2[l] += frac * (v[l] - normscal * vp[l]) / normvp;
		F1[l] -= frac * ((vp[l] - normscal * v[l]) / normv + (v[l] - normscal * vp[l]) / normvp);
	}
}

static inline void updateForces(double aux, const double Kappa, const double factor, double (*const F)[3],
								const double (*const r)[3], const int i, const int j)
{
	aux = exp(-Kappa * aux) * (Kappa + 1. / aux) / (aux * aux) * factor; // 1 g/ps^2/mol
	for (int l = 0; l < 3; l++)
	{
		const double diff = (r[i][l] - r[j][l]) * aux; // 1 g*A/ps^2/mol
		#pragma omp atomic
		F[i][l] += diff;
		#pragma omp atomic
		F[j][l] -= diff;
	}
}

void DH(double (*F)[3], double (*r)[3], double q, double IonicStrength, double Eps, int n)
{
	/*
	Debye-Hückel equation Forces between n particles
	eps = 80.0;
	kappa = 0.103; 			// unit: 1/A	taken from Arya 2014 paper (usually 0.033 for 10mM salt) Mueller 0.103 for
	100mM salt kappa^2 proportional to I f_pi_eps0 = 1.1126*pow(10,-10);

	dh = q1 * 1.6 * pow(10, -19) * q2 * 1.6 * pow(10, -19) / (f_pi_eps0 * eps) * 1 / r * pow(10, 10) * exp(-kappa * r) /
	(1.6 * pow(10, -19));	// 10^10 wegen r in Angstrom
	*/
	const double q1 = q, q2 = q;

	const double Eps0 = 1.1126e-10; /*C^2*s^2/(kg*m^3)*/
	const double NA = 6.02214076e23;
	const double Kappa100 = 0.103;
	const double Kappa = Kappa100 * sqrt(IonicStrength / 100);
	const int nn = NN;

	const double factor = q1 * q2 * 1 / (Eps0 * Eps) * 1e9 * NA;

	/*NEW VERSION WITHOUT CRITICAL*/
	/*Separate distance algorithm INTRA*/
	#pragma omp for schedule(static) label("forces-DH-1")
	for (int diagonal = 1 + nn; diagonal < n / 2 - nn; diagonal++)
		// #pragma omp single label("forces-DH-1")
		for (int index = 0; index < n / 2; index++)
		{
			int j = index + diagonal;
			int i = index;

			if (j >= n / 2)
			{
				i = j;
				j = index + n / 2;
			}

			const double sqDist = distSq(r[i], r[j]);
			if (sqDist < 30*30)
				updateForces(sqrt(sqDist), Kappa, factor, F, r, i, j);
		}

	/*Separate distance algorithm INTER*/
	#pragma omp for schedule(static) label("forces-DH-2")
	for (int diagonal = 1 + nn; diagonal < n / 2 - nn; diagonal++)
	{
		// #pragma omp single label("forces-DH-2")
		for (int i = 0; i < n / 2; i++)
		{
			int j = n - diagonal - i - 1;

			if (j < n / 2)
				j = j + n / 2;

			const double sqDist = distSq(r[i], r[j]);
			if (sqDist < 30*30)
				updateForces(sqrt(sqDist), Kappa, factor, F, r, i, j);
		}
	}
}

void LJ(double (*F)[3], double (*r)[3], double EpsLJ, double sigma, int n)
{
	/*Lennard-Jones equation, with distance sigma, and constant EpsLJ*/
	int i, j, l;
	double aux, aux2, ratio;

	for (j = 0; j < n; j++)
		for (i = j + 1; i < n; i++)
		{
			const double sqDist = distSq(r[i], r[j]);
			ratio = sigma / sqrt(sqDist);
			ratio = ratio * ratio;		   /*1+1=2*/
			ratio = ratio * ratio * ratio; /*2+2+2=6*/
			aux = 24 * EpsLJ * (2 * ratio * ratio - ratio) / sqDist;
			for (l = 0; l < 3; l++)
			{
				aux2 = (r[i][l] - r[j][l]) * aux;
				F[i][l] += aux2;
				F[j][l] -= aux2;
			}
		}
}

void Brownian(double sigma, double *mass, double (*F)[3], int n)
{
	#pragma omp for label("forces-brownian")
	for (int i = 0; i < n; i++)
	{
		for (int l = 0; l < 3; l++)
		{
			F[i][l] += sigma / sqrt(mass[i]) * sampleNormal();
		}
	}
}

double BondEnergyij(double (*r)[3], double d, double *K, int i, int j)
{
	/*Applies BondEnergy to the corresponding Points*/
	return BondEnergy(r[i], r[j], d, K);
}

double BondEnergy(double r[3], double s[3], double d, double K[Order])
{
	/*Computes the Energy associated with a Bonded interaction between 2 particles,
	with equilibrium distance d and polynomial Constants K*/
	int j;
	double sum = 0., aux, aux2;
	aux = dist(r, s);
	for (j = Order, aux2 = 0.; j > 0; j--)
		aux2 = (aux2 + K[j - 1]) * (aux - d);
	sum += aux2 * (aux - d);
	return sum;
}

double AngleEnergyijk(double (*r)[3], double angle, double *K, int i, int j, int k)
{
	/*Applies AngleEnergy to the corresponding Points*/
	return AngleEnergy(r[i], r[j], r[k], angle, K);
}

double AngleEnergy(double r0[3], double r1[3], double r2[3], double angle, double K[Order])
{
	/*Computes the Energy associated with an Angle interaction between 3 particles,
	with equilibrium angle 'angle' and polynomial Constants K*/
	int j;
	double sum = 0., aux, aux2;

	aux = PointsToAngle(r0, r1, r2);
	for (j = Order, aux2 = 0.; j > 0; j--)
		aux2 = (aux2 + K[j - 1]) * (aux - angle);
	sum += aux2 * (aux - angle);
	return sum;
}

void DHEnergy(double (*r)[3], double q, double IonicStrength, double Eps, int n, double *E)
{
	/*
	Debye-Hückel equation Energies
	eps = 80.0;
	kappa = 0.103; 			// unit: 1/A	taken from Arya 2014 paper (usually 0.033 for 10mM salt) Mueller 0.103 for
	100mM salt kappa^2 proportional to I f_pi_eps0 = 1.1126*pow(10,-10);

	dh = q1 * 1.6 * pow(10, -19) * q2 * 1.6 * pow(10, -19) / (f_pi_eps0 * eps) * 1 / r * pow(10, 10) * exp(-kappa * r) /
	(1.6 * pow(10, -19));	// 10^10 wegen r in Angstrom
	*/
	const double q1 = q, q2 = q;

	const double Eps0 = 1.1126e-10; /*C^2*s^2/(kg*m^3)*/
	const double NA = 6.02214076e23;
	const double Kappa100 = 0.103;
	const double Kappa = Kappa100 * sqrt(IonicStrength / 100);
	int nn = NN;

	#pragma omp single label("forces-DHEnergy-zero")
	*E = 0.;

	const double factor = q1 * q2 * 1 / (Eps0 * Eps) * 1e9 * NA;

	#pragma omp for schedule(dynamic) reduction(+:E[0:1]) label("forces-DHEnergy")
	for (int j = 0; j < n; j++)
	{
		double tE = 0;
		for (int i = j + 3; i < n; i++)
		{
			if (i + j <= n - 1 - nn || i + j >= n - 1 + nn)
			{
				const double dSq = distSq(r[i], r[j]);
				if (dSq < 30*30)
				{
					const double aux = sqrt(dSq);
					tE += exp(-Kappa * aux) / aux * factor; /*1 g/ps^2/mol*/
				}
			}
		}
		*E += tE;
	}
}

double LJEnergy(double (*r)[3], double EpsLJ, double sigma, int n)
{
	/*Lennard-Jones equation, with distance sigma, and constant EpsLJ*/
	int i, j;
	double aux, ratio;

	double E = 0.;

	for (j = 0; j < n; j++)
		for (i = j + 1; i < n; i++)
		{
			aux = dist(r[i], r[j]);
			ratio = sigma / aux;
			ratio = ratio * ratio;		   /*1+1=2*/
			ratio = ratio * ratio * ratio; /*2+2+2=6*/
			aux = 24 * EpsLJ * (2 * ratio * ratio - ratio) / (aux * aux);
			E += 4 * EpsLJ * (ratio * ratio - ratio);
		}
	return E;
}

double kineticEnergy(double (*v)[3], int n, double *m)
{
	/*Computes the kinetic energy of a system with n particles, velocities v, and masses m*/
	int i;
	double sum = 0.;
	for (i = 0; i < n; i++)
		sum += m[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

	return sum / 2.;
}
