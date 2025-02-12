#include "settings.h"

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

double distSq(double *p, double *q);
double dist(double *p, double *q);

double norm(double *p);
double scalar(double *u, double *v);
void prodVect(double *u, double *v, double *w);

double PointsToAngle(double *r0, double *r1, double *r2);

/*These 4 functions are currently not being used, they were originally design to apply isometries to static objects*/
int FindIsometry(double **pin, double **pfin, double **M);
void ApplyIsometry(double **M, double (*u)[3], double (*v)[3], int n);
void StaticNucleosome(double **p, double **Val, int n);
void RebuildNucleosome(double **q, double **Val, int n);

void seedThreads(const long baseSeed);
double sampleNormal();

#endif // LINEAR_ALGEBRA_H
