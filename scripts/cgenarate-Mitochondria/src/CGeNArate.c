#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#include "TetramersCircular.h"
#include "TetramersLinear.h"
#include "forces.h"
#include "linear_algebra.h"
#include "memory.h"
#include "settings.h"
#include "config.h"

#define SCANF_CHK1(fn)                                                                                                 \
	do                                                                                                                 \
	{                                                                                                                  \
		if (fn != 1)                                                                                                   \
		{                                                                                                              \
			fprintf(stderr, "Error scanf: %s:%d\n", __FILE__, __LINE__);                                               \
			exit(1);                                                                                                   \
		}                                                                                                              \
	} while (0)

double CGsimulation(double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
				 double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[NNLong],
				 double KLong[NNLong][Order], int nn);

double gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec * 1e-6;
}

void usage(const char* name) {
	fprintf(stderr, "Usage: %s <settings.toml>\n", name);
}

int main(int argc, char *argv[])
{

	if (argc < 2) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	config_parse(argv[1]);

	int i, j, k;
	double DTetramers[AmountOfTetramers][ConstantsPerTetramer];		   /*#Tetramers, #Springs/Tetramer*/
	double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order]; /*#Tetramers, #Springs/Tetramer*/
	double DLong[NNLong];											   /*#LongSprings*/
	double KLong[NNLong][Order];									   /*#LongSprings*/

	FILE *Dfile, *Kfile, *DLongfile, *KLongfile;

	Dfile = fopen_read("input/MyD.txt", "Error opening Parameter file\n");
	for (i = 0; i < AmountOfTetramers; i++)
		for (j = 0; j < ConstantsPerTetramer; j++)
			SCANF_CHK1(fscanf(Dfile, "%lf", &DTetramers[i][j]));
	fclose(Dfile);

	Kfile = fopen_read("input/MyK.txt", "Error opening Parameter file\n");
	for (i = 0; i < AmountOfTetramers; i++)
		for (j = 0; j < ConstantsPerTetramer; j++)
			for (k = 0; k < Order; k++)
				SCANF_CHK1(fscanf(Kfile, "%lf", &KTetramers[i][j][k]));
	fclose(Kfile);

	DLongfile = fopen_read("input/MyDLong.txt", "Error opening Parameter file\n");
	for (j = 0; j < NNLong; j++)
		SCANF_CHK1(fscanf(DLongfile, "%lf", &DLong[j]));
	fclose(DLongfile);

	KLongfile = fopen_read("input/MyKLong.txt", "Error opening Parameter file\n");
	for (j = 0; j < NNLong; j++)
		for (k = 0; k < Order; k++)
			SCANF_CHK1(fscanf(KLongfile, "%lf", &KLong[j][k]));
	fclose(KLongfile);

	double took = CGsimulation(DTetramers, KTetramers, DLong, KLong, NN);
	printf("Simulation time: %lf\n", took);

	return 0;
}

double CGsimulation(double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
				 double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[NNLong],
				 double KLong[NNLong][Order], int nn)
{
	int loadrestart = 0, PrintVectors = 0, PlotPrint = 0, Linear = 0;
	double T0 = 0.0, dt0 = 0.0;
	int N = 1, nframes = 0, step = 0;
	char S[200];

	FILE *infile, *energyplot, *plotfile, *coordinatesfile, *rstfile;

	loadrestart = cgen_config.load;
	PrintVectors = cgen_config.printVectors;
	PlotPrint = cgen_config.printPlot;
	Linear = cgen_config.isLinear;
	T0 = cgen_config.T0;
	dt0 = cgen_config.dt0;
	N = cgen_config.N;
	nframes = cgen_config.frames;
	step = cgen_config.steps;

	fprintf(stderr, "SETTINGS:\n");
	fprintf(stderr, "loadrestart = %d\n", loadrestart);
	fprintf(stderr, "PrinttVectors = %d\n", PrintVectors);
	fprintf(stderr, "PlotPrint = %d\n", PlotPrint);
	fprintf(stderr, "Linear = %d\n", Linear);
	fprintf(stderr, "T0 = %lf\n", T0);
	fprintf(stderr, "dt0 = %lf\n", dt0);
	fprintf(stderr, "N = %d\n", N);
	fprintf(stderr, "nframes = %d\n", nframes);
	fprintf(stderr, "step = %d\n", step);

	config_dump(stderr);

	void (*ToTetramers)(char *RNAME, int n, int *TNums);
	void (*MassFromResidue)(char *RNAME, int n, double *mass);
	void (*LinkerForcesTetramers)(double(*F)[3], double(*r)[3], int n, int nn, int *TNums,
								  double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
								  double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4],
								  double KLong[4][Order]);
	void (*EnergyTetramers)(
		double(*r)[3], int n, int nn, int *TNums, double DTetramers[AmountOfTetramers][ConstantsPerTetramer],
		double KTetramers[AmountOfTetramers][ConstantsPerTetramer][Order], double DLong[4], double KLong[4][Order],
		double *energyTet);
	void (*PrintPlots)(double(*r)[3], int n, int nn, FILE *plotfile);
	if (Linear)
	{
		ToTetramers = ToTetramersLinear;
		MassFromResidue = MassFromResidueLinear;
		LinkerForcesTetramers = LinkerForcesTetramersLinear;
		EnergyTetramers = EnergyTetramersLinear;
		PrintPlots = PrintPlotsLinear;
	}
	else
	{
		ToTetramers = ToTetramersCircular;
		MassFromResidue = MassFromResidueCircular;
		LinkerForcesTetramers = LinkerForcesTetramersCircular;
		EnergyTetramers = EnergyTetramersCircular;
		PrintPlots = PrintPlotsCircular;
	}

	/*Indices*/
	int i, l;
	const int n = N, NucNum = 0;

	if (NucNum > 0)
	{
		printf("Multiple linkers not implemented\n");
		exit(1);
	}

	/*Constants definition 1 g*A^2/ps^2 = 10 J*/
	const double dt = dt0; /*ps*/
	const double IonicStrength = 100., Eps = 80;
	const double q = 1.6e-19;
	const double gamma = .5; /*1/ps*/
	const double kB = 0.8314463 /* g*A^2/(K*mol*ps^2)=0.001986*kcal/(mol*K)*/;
	const double T = T0 /*K*/;
	const double epsilon = exp(-gamma * dt);
	const double sigma = sqrt(kB * T * (1 - exp(-2 * gamma * dt))) /* sqrt(g/mol)*A/ps^2)*/;

	/*Coordinates data*/
	double (* const point)[N][3] = malloc(sizeof(double[NucNum + 1][N][3]));
	double (* const v)[N][3] = malloc(sizeof(double[NucNum + 1][N][3]));
	double (* const F)[N][3] = malloc(sizeof(double[NucNum + 1][N][3]));

	double(* const mass)[N] = calloc(NucNum + 1, sizeof(double[N])); /*g/mol*/
	int(* const TNum)[N] = calloc(NucNum + 1, sizeof(int[N / 2])); /*Linear DNA uses one less tetramer, which is left empty */
	char(* const RNAME)[N] = calloc(NucNum + 1, sizeof(char[N]));	/*Residue*/

	/*Auxiliary and temporary arrays*/
	char aname[4 + 1], res[10];
	char auxstr[8 + 1];

	time_t mytime = time(NULL);
	mytime = 24;
	printf("TEST RUN ");
	printf("%d\n", (int)mytime);
	seedThreads(mytime);

	energyplot = cgen_config.energyFile;
	coordinatesfile = cgen_config.outputFile;

	plotfile = cgen_config.plotsFile;

	/*Read structure*/
	infile = cgen_config.inputFile;

	i = 0;
	while (fgets(S, sizeof(S), infile) != NULL)
	{
		/*Read Atom Name. Is it C1'?*/
		sprintf(aname, "%.*s", 4, S + 12);

		if (!strcmp(aname, " C1'"))
		{
			/*Read Coordinates*/
			for (l = 0; l < 3; l++)
			{
				sprintf(auxstr, "%.*s", 8, S + 30 + 8 * l);
				point[0][i][l] = strtod(auxstr, NULL);
			}

			/*Read Residue Name*/
			sprintf(res, "%.*s", 3, S + 17);
			if (res[0] == ' ')
				RNAME[0][i] = res[2];
			else
				RNAME[0][i] = res[1];
			i++;
		}
		memset(S, 0, strlen(S));
	}
	fclose(infile);

	ToTetramers(RNAME[0], N / 2, TNum[0]);
	MassFromResidue(RNAME[0], N, mass[0]);

	fprintf(coordinatesfile, "Cpptraj Generated trajectory                                                    \n");
	PrintCoordinates(point, n, coordinatesfile);

	/*Iterate the system*/
	if (loadrestart == 2)
	{
		rstfile = cgen_config.restartFile;

		i = 0;
		while (fgets(S, sizeof(S), rstfile) != NULL)
		{
			for (l = 0; l < 3; l++)
			{
				sprintf(auxstr, "%.*s", 8, S + 30 + 8 * l);
				v[0][i][l] = strtod(auxstr, NULL);
			}
			i++;
		}
		fclose(rstfile);
	}

	double vtot[3];
	double energyDH, energyTet, energyV;

	double t0 = gettime();

	const double  sqKBT = sqrt(kB * T);

	#ifdef DEBUG
	double *frameTimes = calloc(nframes, sizeof(double));
	#else
	double *frameTimes = NULL;
	#endif


	#pragma omp parallel default(none) shared(point, v, F, TNum, mass, RNAME) shared(DTetramers, KTetramers, DLong, KLong) \
	shared(n, nn, NucNum, dt, step, nframes) shared(Eps, IonicStrength, q, epsilon, sigma, sqKBT)                      \
	shared(coordinatesfile, energyplot, loadrestart, PlotPrint, plotfile)                        \
	shared(ToTetramers, MassFromResidue, LinkerForcesTetramers, EnergyTetramers, PrintPlots) shared(vtot)              \
	shared(energyDH, energyTet, energyV)										\
	shared(frameTimes)
	{
		if (loadrestart != 2)
		{
			#pragma omp single label("CGeNArate-brownian-pre")
			memset(v, 0, sizeof(double[NucNum + 1][n][3]));
			for (int r = 0; r < NucNum + 1; r++)
			{
				Brownian(sqKBT, mass[r], v[r], n);
			}
		}

		DHEnergy(point[0], q, IonicStrength, Eps, n, &energyDH);
		EnergyTetramers(point[0], n, nn, TNum[0], DTetramers, KTetramers, DLong, KLong, &energyTet);

		#pragma omp single label("CGeNArate-print-energy-init")
		{
			const double kinetic_init = kineticEnergy(v[0], n, mass[0]);
			fprintf(energyplot, "%lf,%lf,%lf\n", energyDH, energyTet, kinetic_init);
		}

		/*Starting force*/
		for (int r = 0; r < NucNum + 1; r++)
		{
			LinkerForcesTetramers(F[r], point[r], n, nn, TNum[r], DTetramers, KTetramers, DLong, KLong);
			DH(F[r], point[r], q, IonicStrength, Eps, n);
		}

		#pragma omp single label("CGeNArate-print-plots-init")
		if (PlotPrint)
			PrintPlots(point[0], n, nn, plotfile);

	/*Full process*/
		for (int k = 0; k < nframes; k++)
		{
			#ifdef DEBUG
			#pragma omp single label("CGeNArate-frame-debug-pre")
			{
				frameTimes[k] = gettime();
			}
			#endif

			for (int h = 0; h < step; h++)
			{

				#ifdef DEBUG
				#pragma omp single label("CGeNArate-frame-debug-step")
				{
					fprintf(stderr, "\rFrame %d : %d/%d       ", k, h, step);
				}
				#endif

				#pragma omp single label("CGeNArate-for-step-1")
				for (int r = 0; r < NucNum + 1; r++)
					for (int i = 0; i < n; i++)
						// printf("%d %d\n", BlockNum, i);
						for (int l = 0; l < 3; l++)
						{
							v[r][i][l] += dt / 2. * F[r][i][l] / mass[r][i];
							point[r][i][l] += dt * v[r][i][l];
						}

				// #pragma omp single //IT IS INSIDE THE FUNCTION
				for (int r = 0; r < NucNum + 1; r++)
				{
					LinkerForcesTetramers(F[r], point[r], n, nn, TNum[r], DTetramers, KTetramers, DLong, KLong);
					DH(F[r], point[r], q, IonicStrength, Eps, n);
				}

				#pragma omp single label("CGeNArate-for-step-2")
				for (int r = 0; r < NucNum + 1; r++)
					for (int i = 0; i < n; i++)
						for (int l = 0; l < 3; l++)
						{
							v[r][i][l] += dt / 2. * F[r][i][l] / mass[r][i];
							v[r][i][l] = epsilon * v[r][i][l];
						}

				// #pragma omp single //IT IS INSIDE THE FUNCTION
				for (int r = 0; r < NucNum + 1; r++)
					Brownian(sigma, mass[r], v[r], n); // Can we do parallel RNG?

				#pragma omp single label("CGeNArate-step-end")
				{
					memset(vtot, 0, sizeof(vtot));

					for (int r = 0; r < NucNum + 1; r++)
						for (int i = 0; i < n; i++)
							for (int l = 0; l < 3; l++)
								vtot[l] += v[r][i][l];

					for (int r = 0; r < NucNum + 1; r++)
						for (int i = 0; i < n; i++)
							for (int l = 0; l < 3; l++)
								v[r][i][l] -= vtot[l] / n;
				}
			}

			/*End of Cycle Printing*/
			{
				DHEnergy(point[0], q, IonicStrength, Eps, n, &energyDH);
				EnergyTetramers(point[0], n, nn, TNum[0], DTetramers, KTetramers, DLong, KLong, &energyTet);
				#pragma omp single label("CGeNArate-kineticEnergy")
				energyV = kineticEnergy(v[0], n, mass[0]);
			}

			#ifdef DEBUG
			#pragma omp single label("CGeNArate-frame-debug-end")
			{
				frameTimes[k] = gettime() - frameTimes[k];
				fprintf(stderr, "\rFrame %d : %d/%d DONE in %lf\n", k, step, step, frameTimes[k]);
			}
			#endif

			#pragma omp sections
			{
				#pragma omp section
				if (PlotPrint)
				{
					PrintPlots(point[0], n, nn, plotfile);
				}

				#pragma omp section
				{
					fprintf(energyplot, "%lf,%lf,%lf\n", energyDH, energyTet, energyV);
				}

				#pragma omp section
				{
					PrintCoordinates(point, n, coordinatesfile);
				}

				#pragma omp section
				if (loadrestart)
				{
					if (freopen(NULL, "w", cgen_config.restartFile) == NULL) {
						fprintf(stderr, "ERROR reopening file: %s\n", strerror(errno));
					}
					for (int i = 0; i < n; i++)
						fprintf(cgen_config.restartFile, "ATOM  %5d %4s %3s  %4d    %8.3lf%8.3lf%8.3lf\n", i + 1, "N ", "A", i + 1,
								v[0][i][0], v[0][i][1], v[0][i][2]);
					fprintf(cgen_config.restartFile, "END\n");
					fflush(cgen_config.restartFile);
				}
				#pragma omp section
				if (loadrestart)
				{
					if (freopen(NULL, "w", cgen_config.currentFile) == NULL) {
						fprintf(stderr, "ERROR reopening file: %s\n", strerror(errno));
					}
					for (int i = 0; i < n; i++)
						fprintf(cgen_config.currentFile, "ATOM  %5d %4s %3c  %4d    %8.3lf%8.3lf%8.3lf\n", i + 1, " C1'", RNAME[0][i],
								i + 1, point[0][i][0], point[0][i][1], point[0][i][2]);
					fprintf(cgen_config.currentFile, "END\n");
					fflush(cgen_config.currentFile);
					
				}
			}
		}
	}
	double t1 = gettime();
	printf("\n");

	free(point);
	free(v);
	free(F);
	free(mass);
	free(TNum);
	free(RNAME);

	if (PlotPrint)
		fclose(plotfile);
	fclose(coordinatesfile);
	fclose(energyplot);

	return t1-t0;
}
