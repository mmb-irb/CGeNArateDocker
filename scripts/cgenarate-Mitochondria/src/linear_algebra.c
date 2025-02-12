#include "linear_algebra.h"
#define _DEFAULT_SOURCE
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


double distSq(double *p, double *q)
{
	/*Distance between two points in R3*/
	return (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
}

inline double dist(double *p, double *q)
{
	return sqrt(distSq(p, q));
}

double norm(double *p)
{
	/*Euclidian norm of vector in R3*/
	return sqrt((p[0]) * (p[0]) + (p[1]) * (p[1]) + (p[2]) * (p[2]));
}

double scalar(double *u, double *v)
{
	/*Scalar product between two vectors in R3*/
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

void prodVect(double *u, double *v, double *w)
{
	/*Vector product between two vectors in R3*/
	w[0] = u[1] * v[2] - u[2] * v[1];
	w[1] = u[2] * v[0] - u[0] * v[2];
	w[2] = u[0] * v[1] - u[1] * v[0];
}

double PointsToAngle(double *r0, double *r1, double *r2)
{
	/*Angle between three points in R3*/
	double v[3], vp[3];
	int l;

	for (l = 0; l < 3; l++)
	{
		v[l] = r0[l] - r1[l];
		vp[l] = r2[l] - r1[l];
	}
	return acos(scalar(v, vp) / (norm(v) * norm(vp)));
}

/*void BrownianSingle(double sigma, double (*F)[3], int n){
	double w = 2.;
	double rx[3];

	int i,l;

	for(i=0; i<n; i++){
		w = 2.;
		while ( w >= 1 ){
			for(l=0; l<3; l++)
				rx[l] = (rand()*2./(RAND_MAX)) - 1.;

			w = rx[0]*rx[0]+rx[1]*rx[1]+rx[2]*rx[2];
		}
		w = sqrt(-2.*log(w)/w);

		for(l=0; l<3; l++){
			F[i][l]+=sigma*w*rx[l];
		}
	}
}*/

static struct drand48_data *g_drandBuff;

// Seeds each thread with <baseSeed> + <THREAD ID>
void seedThreads(const long baseSeed)
{
	#ifdef _OPENMP
	const int MAX_THREADS = omp_get_max_threads();
	#else
	const int MAX_THREADS = 1;
	#endif
	g_drandBuff = malloc(sizeof(struct drand48_data) * MAX_THREADS);

	fprintf(stderr, "Seeding %d threads\n", MAX_THREADS);
	for (int i = 0; i < MAX_THREADS; ++i)
	{
		srand48_r(baseSeed + i, &g_drandBuff[i]);
	}
}

// WARN: seedThreads has to have been called before.
double sampleNormal()
{
	/*random sample of Normal(0,1) following Box-Muller transform*/
	#ifdef _OPENMP
	const int tid = omp_get_thread_num();
	#else
	const int tid = 0;
	#endif

	double u, v;
	drand48_r(&g_drandBuff[tid], &u);
	drand48_r(&g_drandBuff[tid], &v);
	u = u * 2 - 1;
	v = v * 2 - 1;
	const double r = u * u + v * v;
	if (r == 0 || r > 1)
		return sampleNormal();
	const double c = sqrt(-2 * log(r) / r);
	return u * c;
}

int FindIsometry(double **pin, double **pfin, double **M)
{
	/*In R3, Find the 4x4 matrix that transforms the 4 points from pin into the 4 points from pfin*/
	int i, j, k, r;
	double **pinv, *aux;
	double mul;

	pinv = malloc(4 * sizeof(double *));
	if (pinv == NULL)
	{
		printf("Error allocating memory\n");
		exit(1);
	}
	for (i = 0; i < 4; i++)
	{
		pinv[i] = malloc(4 * sizeof(double));
		if (pinv[i] == NULL)
		{
			printf("Error allocating memory\n");
			exit(1);
		}
	}

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			pinv[i][j] = 0;

	for (i = 0; i < 4; i++)
		pinv[i][i] = 1;

	for (k = 0; k < 4; k++)
	{
		for (i = k + 1, r = k; i < 4; i++)
			if (fabs(pin[i][k]) > fabs(pin[r][k]))
				r = i;
		if (r != k)
		{
			aux = pin[r];
			pin[r] = pin[k];
			pin[k] = aux;
			aux = pinv[r];
			pinv[r] = pinv[k];
			pinv[k] = aux;
		}
		if (fabs(pin[k][k]) < tol)
		{
			printf("Error!!!!\n");
		}

		mul = pin[k][k];
		for (j = 0; j < 4; j++)
			pinv[k][j] /= mul;
		for (j = 0; j < 4; j++)
			pin[k][j] /= mul;

		for (i = k + 1; i < 4; i++)
		{
			mul = pin[i][k];
			for (j = 0; j < 4; j++)
				pin[i][j] -= mul * pin[k][j];
			for (j = 0; j < 4; j++)
				pinv[i][j] -= mul * pinv[k][j];
		}
	}
	for (k = 3; k > 0; k--)
	{
		for (i = k - 1; i >= 0; i--)
		{
			mul = pin[i][k];
			for (j = 0; j < 4; j++)
				pin[i][j] -= mul * pin[k][j];
			for (j = 0; j < 4; j++)
				pinv[i][j] -= mul * pinv[k][j];
		}
	}

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			M[i][j] = 0;
			for (k = 0; k < 4; k++)
				M[i][j] += pfin[i][k] * pinv[k][j];
		}

	for (i = 0; i < 4; i++)
		free(pinv[i]);
	free(pinv);

	return 0;
}

void ApplyIsometry(double **M, double (*u)[3], double (*v)[3], int n)
{
	/*Transforms n 3dvectors from u, into n 3dvectors stored in v, according to the isometry M,
	Can be used with u==v, i.e. overwritting vectors in u*/
	int i, j;
	double mov[3];
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < 3; i++)
			mov[i] = M[i][0] * u[j][0] + M[i][1] * u[j][1] + M[i][2] * u[j][2] + M[i][3];
		for (i = 0; i < 3; i++)
			v[j][i] = mov[i];
	}
}

void StaticNucleosome(double **p, double **Val, int n)
{
	/*Computes the relationship between the 4 3dpoints describing a nucleosome, and stores them in the Val vector,
	Can be used with multiple nucleosomes at once, with n>1, and p being an nx3 vector*/
	int i, l;
	double u[3], v[3], w[3];

	for (i = 0; i < n; i++)
	{
		for (l = 0; l < 3; l++)
			u[l] = p[3 + i * 4][l] - p[0 + i * 4][l]; /*v14*/
		Val[i][0] = sqrt(scalar(u, u));
		if (Val[i][0] < tol)
			printf("Plane\n");
		for (l = 0; l < 3; l++)
			u[l] = u[l] / Val[i][0]; /*v14 normalized*/

		for (l = 0; l < 3; l++)
			v[l] = p[2 + i * 4][l] - p[0 + i * 4][l]; /*v13*/
		Val[i][1] = scalar(u, v);
		for (l = 0; l < 3; l++)
			v[l] -= Val[i][1] * u[l]; /*v13 projection over v14Ort*/
		Val[i][2] = sqrt(scalar(v, v));
		if (Val[i][2] < tol)
			printf("Plane\n");
		for (l = 0; l < 3; l++)
			v[l] = v[l] / Val[i][2]; /*v14Ort normalized*/

		for (l = 0; l < 3; l++)
			w[l] = p[1 + i * 4][l] - p[0 + i * 4][l]; /*v12*/
		Val[i][3] = scalar(w, u);
		Val[i][4] = scalar(w, v);

		for (l = 0; l < 3; l++)
			w[l] -= Val[i][3] * u[l] + Val[i][4] * v[l]; /*v12 projection over uxv*/
		Val[i][5] = sqrt(scalar(w, w));					 /*Sign?*/
		if (Val[i][5] < tol)
			printf("Plane\n");
	}
}

void RebuildNucleosome(double **p, double **Val, int n)
{
	/*Rebuilds a nucleosome according to the Val vector, and the 4 3dpoints describing the nucleosome
	Can be used with multiple nucleosomes at once, with n>1, and p being an nx3 vector
	Uses Gramm-Schmidt process*/
	int i, j, l;
	double aux;
	double u[3], v[3], w[3];
	double center[3], offset[3];

	for (i = 0; i < n; i++)
	{

		for (l = 0; l < 3; l++)
		{
			center[l] = 0;
			for (j = 0; j < 4; j++)
				center[l] += p[j + i * 4][l];
		}

		for (l = 0; l < 3; l++)
			u[l] = p[3 + i * 4][l] - p[0 + i * 4][l]; /*v14*/
		aux = sqrt(scalar(u, u));
		if (aux < sqrt(tol))
			printf("Plane\n");
		for (l = 0; l < 3; l++)
			u[l] = u[l] / aux; /*v14 normalized*/
		for (l = 0; l < 3; l++)
			p[3 + i * 4][l] = p[0 + i * 4][l] + Val[i][0] * u[l]; /*p4*/

		for (l = 0; l < 3; l++)
			v[l] = p[2 + i * 4][l] - p[0 + i * 4][l]; /*v13*/
		aux = scalar(u, v);
		for (l = 0; l < 3; l++)
			v[l] -= aux * u[l]; /*v13 projection over v14Ort*/
		aux = sqrt(scalar(v, v));
		if (aux < sqrt(tol))
			printf("Plane\n");
		for (l = 0; l < 3; l++)
			v[l] = v[l] / aux; /*v14Ort normalized*/
		for (l = 0; l < 3; l++)
			p[2 + i * 4][l] = p[0 + i * 4][l] + Val[i][1] * u[l] + Val[i][2] * v[l]; /*p3*/

		prodVect(u, v, w);

		for (l = 0; l < 3; l++)
			p[1 + i * 4][l] = p[0 + i * 4][l] + Val[i][3] * u[l] + Val[i][4] * v[l] + Val[i][5] * w[l]; /*p2*/

		for (l = 0; l < 3; l++)
		{
			offset[l] = 0;
			for (j = 0; j < 4; j++)
				offset[l] += p[j + i * 4][l];
			offset[l] = (center[l] - offset[l]) / 4.;
		}
		for (j = 0; j < 4; j++)
			for (l = 0; l < 3; l++)
				p[j + i * 4][l] += offset[l];
	}
}
