#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <toml.h>

#ifndef CONFIG_H
#define CONFIG_H

#define __internal __attribute__((visibility("hidden")))

enum loadRestart
{
	NOTHING,
	PRINT,
	LOAD_AND_PRINT
};

typedef struct cgen_config
{
	// Restore state from file:
	// - NOTHING -> Does nothing
	// - PRINT -> Only print, does not load
	// - LOAD_AND_PRINT -> Loads and prints
	// Requires restoreFile to be set
	enum loadRestart load;

	// Prints al vector velocities
	bool printVectors;
	// Print distances and angles relevant to the (bonded) Hamiltonian
	bool printPlot;

	// Number of beads
	unsigned int N;

	double T0; // Initial temperature

	double dt0; // TimeStep in picoseconds

	// We record a frame every STEP*dt0 picoseconds
	unsigned int frames;
	unsigned int steps; // Number of steps of duration dt0 on a frame

	bool isLinear;

	// Seed used for RNG.
	unsigned long seed;

	// Input files
	FILE *restartFile;
	FILE *currentFile;

	FILE *inputFile; // pdb file

	// Outputs
	FILE *outputFile; // coordinates
	FILE *energyFile;
	FILE *plotsFile;
} cgen_config_t;

static inline void config_init();
void config_parse(const char *const config_path);
void config_process_env();

void config_restore(toml_table_t *);
void config_output(toml_table_t *);
void config_input(toml_table_t *);
void config_simulation(toml_table_t *);

void config_dump(FILE *fd);

__internal extern cgen_config_t cgen_config;

#endif // CONFIG_H
