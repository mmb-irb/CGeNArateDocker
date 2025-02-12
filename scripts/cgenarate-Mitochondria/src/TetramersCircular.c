#include "TetramersCircular.h"
#include "forces.h"
#include "linear_algebra.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ToTetramersCircular(char *RNAME, int nbp, int *TNums)
{
	/*Reads the bases RNAME associated with each particle, and returns the Tetramer identifier TNums associated with
	 each grup of 4 nucleotides. nbp is the number of basePairs. Tetramers wrap around the circle. Note: this code reads
	 ONLY the Watson strand, and therefore assumes canonical basepairs A-T, C-G*/
	int i, k;
	int aux;
	int num;

	/*Beggining*/
	i = 0;
	num = 0;
	k = nbp - 1;
	aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
	num = num << 2;
	num = num ^ aux;
	for (k = 0; k < 3; k++)
	{
		aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
		num = num << 2;
		num = num ^ aux;
	}

	TNums[i] = num;

	/*Full*/
	for (i = 1; i < nbp - 2; i++)
	{
		for (num = 0, k = -1; k < 3; k++)
		{

			aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
			num = num << 2;
			num = num ^ aux;
		}

		TNums[i] = num;
	}

	/*End*/
	i = nbp - 2;
	for (num = 0, k = -1; k < 2; k++)
	{
		aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
		num = num << 2;
		num = num ^ aux;
	}
	k = 2 - nbp;
	aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
	num = num << 2;
	num = num ^ aux;

	TNums[i] = num;

	/*Connecting tetramer*/
	i = nbp - 1;

	for (num = 0, k = -1; k < 1; k++)
	{
		aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
		num = num << 2;
		num = num ^ aux;
	}
	for (k = 1 - nbp; k < 3 - nbp; k++)
	{
		aux = (RNAME[i + k] == 'C') + 2 * (RNAME[i + k] == 'G') + 3 * (RNAME[i + k] == 'T');
		num = num << 2;
		num = num ^ aux;
	}
	TNums[i] = num;

	return;
}

void MassFromResidueCircular(char *RNAME, int n, double *mass)
{
	/*Gives each particle a mass value corresponding to the nucleotide they represent*/
	double masses[4] = {331.222, 307.197, 347.2243, 322.20834}; /*ACGT*/
	int i;
	for (i = 0; i < n; i++)
		mass[i] = masses[0] * (RNAME[i] == 'A') + masses[1] * (RNAME[i] == 'C') + masses[2] * (RNAME[i] == 'G') +
				  masses[3] * (RNAME[i] == 'T');
}

void LinkerForcesTetramersCircular(double (*F)[3], double (*r)[3], int n, int nn, int *TNums,
								   double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
								   double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
								   double KLong[4][Order])
{

	/*Computes all the forces asociated to Bonded potentials of the quasi-harmonic fan model of CGeNArate.
	There are n particles, with positions r and forces F, with n/2 tetramer identifiers
	The dictionaries of Distances and Constants are DTetramers and KTetramers for the 4mer region,
	and Dlong and Klong for the distant region, nn is the maximum bonded interaction, nn=5*/

	/*The - implies Crick strand (Watson otherwise). A implies first appearance (in the Watson direction), B second*/
	/*0 , 1    , 2     , 3     , 4      , 5      , 6 , 7  , 8  , 9 , 10 , 11 , 12, 13 , 14 , 15
	Bond, Bond-, AngleA, AngleB, Angle-A, Angle-B, 3m, 2mA, 2mB, 1m, bpA, bpB, 1p, 2pA, 2pB, 3p*/

	/*    0   , 1    , 3     , 5      , 8  , 9 , 10 , 11 , 12, 14
	+End: Bond, Bond-, AngleB, Angle-B, 2mB, 1m, bpA, bpB, 1p, 2pB
		0   , 1    , 2     , 4      , 7  , 9 , 10 , 11 , 12, 13
	-End: Bond, Bond-, AngleA, Angle-A, 2mA, 1m, bpA, bpB, 1p, 2pA*/
	int i;

	#pragma omp single label("TCirc-linker-init")
	{
		memset(F, 0, sizeof(double) * n * 3); // circular-1

		/*Beginning*/
		i = 0;
		/*Bond, Bond-*/
		BondForceij(F, r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
		BondForceij(F, r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

		/*AngleA, AngleB, Angle-A, Angle-B*/
		AngleForceijk(F, r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1 + (n / 2), i, i + 1);
		AngleForceijk(F, r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2);
		AngleForceijk(F, r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i - (n / 2), n - 1 - i, n - 2 - i);
		AngleForceijk(F, r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i, n - 3 - i);

		/*3m, 2mA, 2mB, 1m*/
		BondForceij(F, r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1 + (n / 2), n - i - 3);
		BondForceij(F, r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1 + (n / 2), n - i - 2);
		BondForceij(F, r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3);
		BondForceij(F, r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

		/*bpA, bpB*/
		BondForceij(F, r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

		/*1p, 2pA, 2pB, 3p*/
		BondForceij(F, r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i - (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2, n - i - (n / 2));
	}

	/*Tetramers*/
	const int skipSize = 4;
	for (int skipConnections = 0; skipConnections < skipSize; skipConnections++)
		#pragma omp for schedule(static) label("TCirc-linker-main")
		for (int i = 1 + skipConnections; i < n / 2 - 2; i += skipSize)
		{
			/*Bond, Bond-*/
			BondForceij(F, r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
			BondForceij(F, r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

			/*AngleA, AngleB, Angle-A, Angle-B*/
			AngleForceijk(F, r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1);
			AngleForceijk(F, r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2);
			AngleForceijk(F, r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i);
			AngleForceijk(F, r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i, n - 3 - i);

			/*3m, 2mA, 2mB, 1m*/
			BondForceij(F, r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3);
			BondForceij(F, r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2);
			BondForceij(F, r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3);
			BondForceij(F, r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

			/*bpA, bpB*/
			BondForceij(F, r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
			BondForceij(F, r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

			/*1p, 2pA, 2pB, 3p*/
			BondForceij(F, r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
			BondForceij(F, r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i);
			BondForceij(F, r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2, n - i - 1);
			BondForceij(F, r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2, n - i);
		}

	#pragma omp single label("TCirc-linker-end")
	{
		/*End*/
		i = n / 2 - 2;
		/*Bond, Bond-*/
		BondForceij(F, r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
		BondForceij(F, r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

		/*AngleA, AngleB, Angle-A, Angle-B*/
		AngleForceijk(F, r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1);
		AngleForceijk(F, r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2 - (n / 2));
		AngleForceijk(F, r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i);
		AngleForceijk(F, r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i,
					  n - 3 - i + (n / 2));

		/*3m, 2mA, 2mB, 1m*/
		BondForceij(F, r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3 + (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2);
		BondForceij(F, r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3 + (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

		/*bpA, bpB*/
		BondForceij(F, r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

		/*1p, 2pA, 2pB, 3p*/
		BondForceij(F, r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i);
		BondForceij(F, r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2 - (n / 2), n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2 - (n / 2), n - i);

		/*Connecting*/
		i = n / 2 - 1;
		/*Bond, Bond-*/
		BondForceij(F, r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1 - (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i + (n / 2));

		/*AngleA, AngleB, Angle-A, Angle-B*/
		AngleForceijk(F, r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1 - (n / 2));
		AngleForceijk(F, r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1 - (n / 2), i + 2 - (n / 2));
		AngleForceijk(F, r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i + (n / 2));
		AngleForceijk(F, r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i + (n / 2),
					  n - 3 - i + (n / 2));

		/*3m, 2mA, 2mB, 1m*/
		BondForceij(F, r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3 + (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2 + (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3 + (n / 2));
		BondForceij(F, r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2 + (n / 2));

		/*bpA, bpB*/
		BondForceij(F, r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1 - (n / 2), n - i - 2 + (n / 2));

		/*1p, 2pA, 2pB, 3p*/
		BondForceij(F, r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1 - (n / 2), n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1 - (n / 2), n - i);
		BondForceij(F, r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2 - (n / 2), n - i - 1);
		BondForceij(F, r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2 - (n / 2), n - i);

		/*LONG*/
		for (int i = 0; i < n / 2 - 5; i++)
			BondForceij(F, r, DLong[0], KLong[0], i, n - i - 6);

		for (int i = 0; i < n / 2 - 4; i++)
			BondForceij(F, r, DLong[1], KLong[1], i, n - i - 5);

		/*5m'*/
		i = n / 2 - 5;
		BondForceij(F, r, DLong[0], KLong[0], i, n - i - 6 + (n / 2));
		
		for (i = n / 2 - 4; i < n / 2; i++)
		{
			/*5m', 4m'*/
			BondForceij(F, r, DLong[0], KLong[0], i, n - i - 6 + (n / 2));
			BondForceij(F, r, DLong[1], KLong[1], i, n - i - 5 + (n / 2));
		}

		for (i = 0; i < 4; i++)
		{
			/*4p', 5p'*/
			BondForceij(F, r, DLong[2], KLong[2], i, n - i + 3 - (n / 2));
			BondForceij(F, r, DLong[3], KLong[3], i, n - i + 4 - (n / 2));
		}

		{
			/*5p'*/
			i = 4;
			BondForceij(F, r, DLong[3], KLong[3], i, n - i + 4 - (n / 2));
		}

		for (int i = 4; i < n / 2; i++)
			BondForceij(F, r, DLong[2], KLong[2], i, n - i + 3);

		for (int i = 5; i < n / 2; i++)
			BondForceij(F, r, DLong[3], KLong[3], i, n - i + 4);
	}
}

void EnergyTetramersCircular(double (*r)[3], int n, int nn, int *TNums /*, int* TShifts*/,
							 double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
							 double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
							 double KLong[4][Order], double *const E)
{
	/*Computes all the Energies asociated to Bonded potentials of the quasi-harmonic fan model of CGeNArate.
	There are n particles, with positions r, with n/2 tetramer identifiers
	The dictionaries of Distances and Constants are DTetramers and KTetramers for the 4mer region,
	and Dlong and Klong for the distant region, nn is the maximum bonded interaction, nn=5*/

	/*The - implies Crick strand (Watson otherwise). A implies first appearance (in the Watson direction), B second*/
	/*0 , 1    , 2     , 3     , 4      , 5      , 6 , 7  , 8  , 9 , 10 , 11 , 12, 13 , 14 , 15
	Bond, Bond-, AngleA, AngleB, Angle-A, Angle-B, 3m, 2mA, 2mB, 1m, bpA, bpB, 1p, 2pA, 2pB, 3p*/

	/*    0   , 1    , 3     , 5      , 8  , 9 , 10 , 11 , 12, 14
	+End: Bond, Bond-, AngleB, Angle-B, 2mB, 1m, bpA, bpB, 1p, 2pB
		  0   , 1    , 2     , 4      , 7  , 9 , 10 , 11 , 12, 13
	-End: Bond, Bond-, AngleA, Angle-A, 2mA, 1m, bpA, bpB, 1p, 2pA*/

	int i;

	/*Beginning*/
	#pragma omp single label("TCirc-energy-init")
	{
		*E = 0;
		double tE = 0;
		i = 0;
		/*Bond, Bond-*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

		/*AngleA, AngleB, Angle-A, Angle-B*/
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1 + (n / 2), i, i + 1);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2);
		tE +=
			AngleEnergyijk(r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i - (n / 2), n - 1 - i, n - 2 - i);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i, n - 3 - i);

		/*3m, 2mA, 2mB, 1m*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1 + (n / 2), n - i - 3);
		tE += BondEnergyij(r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1 + (n / 2), n - i - 2);
		tE += BondEnergyij(r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3);
		tE += BondEnergyij(r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

		/*bpA, bpB*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

		/*1p, 2pA, 2pB, 3p*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i - (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2, n - i - (n / 2));
		*E += tE;
	}

	/*Tetramers*/
	#pragma omp for reduction(+:E[0:1]) label("TCirc-energy-main")
	for (i = 1; i < n / 2 - 2; i++)
	{
		/*Bond, Bond-*/
		double tE = 0.;
		tE += BondEnergyij(r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

		/*AngleA, AngleB, Angle-A, Angle-B*/
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i, n - 3 - i);

		/*3m, 2mA, 2mB, 1m*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3);
		tE += BondEnergyij(r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2);
		tE += BondEnergyij(r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3);
		tE += BondEnergyij(r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

		/*bpA, bpB*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

		/*1p, 2pA, 2pB, 3p*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i);
		tE += BondEnergyij(r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2, n - i);

		*E += tE;
	}
	#pragma omp single label("TCirc-energy-end")
	{
		/*End*/
		double tE = 0;
		i = n / 2 - 2;
		/*Bond, Bond-*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);

		/*AngleA, AngleB, Angle-A, Angle-B*/
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2 - (n / 2));
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i);
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i,
							 n - 3 - i + (n / 2));

		/*3m, 2mA, 2mB, 1m*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3 + (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2);
		tE += BondEnergyij(r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3 + (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);

		/*bpA, bpB*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);

		/*1p, 2pA, 2pB, 3p*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i);
		tE += BondEnergyij(r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2 - (n / 2), n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2 - (n / 2), n - i);

		/*Connecting*/
		i = n / 2 - 1;
		/*Bond, Bond-*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1 - (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i + (n / 2));

		/*AngleA, AngleB, Angle-A, Angle-B*/
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1 - (n / 2));
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1 - (n / 2), i + 2 - (n / 2));
		tE +=
			AngleEnergyijk(r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i + (n / 2));
		tE += AngleEnergyijk(r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i + (n / 2),
							 n - 3 - i + (n / 2));

		/*3m, 2mA, 2mB, 1m*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3 + (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2 + (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3 + (n / 2));
		tE += BondEnergyij(r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2 + (n / 2));

		/*bpA, bpB*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1 - (n / 2), n - i - 2 + (n / 2));

		/*1p, 2pA, 2pB, 3p*/
		tE += BondEnergyij(r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1 - (n / 2), n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1 - (n / 2), n - i);
		tE += BondEnergyij(r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2 - (n / 2), n - i - 1);
		tE += BondEnergyij(r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2 - (n / 2), n - i);

		/*LONG*/
		for (i = 0; i < n / 2 - 5; i++)
		{
			/*5m, 4m*/
			tE += BondEnergyij(r, DLong[0], KLong[0], i, n - i - 6);
			tE += BondEnergyij(r, DLong[1], KLong[1], i, n - i - 5);
		}
		{
			/*5m', 4m*/
			i = n / 2 - 5;
			tE += BondEnergyij(r, DLong[0], KLong[0], i, n - i - 6 + (n / 2));
			tE += BondEnergyij(r, DLong[1], KLong[1], i, n - i - 5);
			for (i = n / 2 - 4; i < n / 2; i++)
			{
				/*5m', 4m'*/
				tE += BondEnergyij(r, DLong[0], KLong[0], i, n - i - 6 + (n / 2));
				tE += BondEnergyij(r, DLong[1], KLong[1], i, n - i - 5 + (n / 2));
			}

			for (i = 0; i < 4; i++)
			{
				/*4p', 5p'*/
				tE += BondEnergyij(r, DLong[2], KLong[2], i, n - i + 3 - (n / 2));
				tE += BondEnergyij(r, DLong[3], KLong[3], i, n - i + 4 - (n / 2));
			}
			/*4p, 5p'*/
			i = 4;
			tE += BondEnergyij(r, DLong[2], KLong[2], i, n - i + 3);
			tE += BondEnergyij(r, DLong[3], KLong[3], i, n - i + 4 - (n / 2));
		}

		for (i = 5; i < n / 2; i++)
		{
			/*4p, 5p*/
			tE += BondEnergyij(r, DLong[2], KLong[2], i, n - i + 3);
			tE += BondEnergyij(r, DLong[3], KLong[3], i, n - i + 4);
		}
		*E += tE;
	}
}

void PrintPlotsCircular(double (*r)[3], int n, int nn, FILE *plotfile)
{

	int i, j, k;

	for (i = 0; i < n / 2; i++)
	{

		if (i + 1 < n / 2)
		{
			fprintf(plotfile, "%lf,", dist(r[i], r[i + 1]));
			fprintf(plotfile, "%lf,", dist(r[n - 1 - i], r[n - 2 - i]));

			if (i > 0)
			{
				fprintf(plotfile, "%lf,", PointsToAngle(r[i - 1], r[i], r[i + 1]));
				fprintf(plotfile, "%lf,", PointsToAngle(r[n - i], r[n - 1 - i], r[n - 2 - i]));
			}
			else if (i == 0)
			{ /*i == 0*/
				fprintf(plotfile, "%lf,", PointsToAngle(r[i - 1 + (n / 2)], r[i], r[i + 1]));
				fprintf(plotfile, "%lf,", PointsToAngle(r[n - i - (n / 2)], r[n - 1 - i], r[n - 2 - i]));
			}
			else
			{
				printf("Wrong indices %d,\n", i);
				exit(1);
			}
		}
		else if (i == n / 2 - 1)
		{ /*i == n/2 -1*/
			fprintf(plotfile, "%lf,", dist(r[i], r[i + 1 - (n / 2)]));
			fprintf(plotfile, "%lf,", dist(r[n - 1 - i], r[n - 2 - i + (n / 2)]));

			fprintf(plotfile, "%lf,", PointsToAngle(r[i - 1], r[i], r[i + 1 - (n / 2)]));
			fprintf(plotfile, "%lf,", PointsToAngle(r[n - i], r[n - 1 - i], r[n - 2 - i + (n / 2)]));
		}
		else
		{
			printf("Wrong indices %d,\n", i);
			exit(1);
		}

		for (k = n - i - 1 - nn - 2, j = 2; j < 2 * nn + 3; j++)
			if (k + j < n && k + j > n / 2 - 1)
			{
				fprintf(plotfile, "%lf,", dist(r[i], r[k + j]));
			}
			else if (k + j >= n)
			{
				fprintf(plotfile, "%lf,", dist(r[i], r[k + j - (n / 2)]));
			}
			else if (k + j <= n / 2 - 1)
			{
				fprintf(plotfile, "%lf,", dist(r[i], r[k + j + (n / 2)]));
			}
			else
			{
				printf("Wrong indices %d, %d, %d\n", k, j, i);
				exit(1);
			}
	}
	fprintf(plotfile, "\n");
	return;
}
