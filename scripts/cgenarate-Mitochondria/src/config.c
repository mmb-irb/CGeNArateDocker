#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <toml.h>

#include "config.h"

__internal cgen_config_t cgen_config;

static void error(const char *msg, const char *msg1)
{
	fprintf(stderr, "ERROR: %s%s\n", msg, msg1 ? msg1 : "");
	exit(1);
}

// Initializes default values
static inline void config_init()
{
	cgen_config.load = false;
	cgen_config.printVectors = false;
	cgen_config.printPlot = false;

	cgen_config.seed = 0;
}

void config_parse(const char *const config_path)
{
	FILE *f = fopen(config_path, "r");
	if (f == NULL)
	{
		error("Failed to open config file for reading", NULL);
	}

	config_init();

	char errbuf[200];
	toml_table_t *conf = toml_parse_file(f, errbuf, sizeof(errbuf));
	fclose(f);

	if (!conf)
	{
		error("cannot parse - ", errbuf);
	}

	config_restore(toml_table_in(conf, "restore"));
	config_simulation(toml_table_in(conf, "simulation"));
	config_output(toml_table_in(conf, "output"));
	config_input(toml_table_in(conf, "input"));

	config_process_env();

	toml_free(conf);
}

static inline char *getFullKey(toml_table_t *t, const char *const key)
{
	const char *tablekey = toml_table_key(t);
	char *fullKey = malloc(strlen(tablekey) + strlen(key) + 2);
	strcpy(fullKey, tablekey);
	strcat(fullKey, ".");
	strcpy(fullKey, key);
	return fullKey;
}

// Reads filename and attempts to open in the given mode.
static inline FILE *cgen_toml_file_in(toml_table_t *t, const char *const key, const char *const mode)
{
	toml_datum_t path = toml_string_in(t, key);
	if (!path.ok)
	{
		error("Key not found", getFullKey(t, key));
	}
	FILE *file = fopen(path.u.s, mode);
	if (file == NULL)
		error("Failed to open file at ", path.u.s);
	free(path.u.s);
	return file;
}

void config_input(toml_table_t *input)
{
	if (input == NULL)
		error("input section missing in TOML", "");

	cgen_config.inputFile = cgen_toml_file_in(input, "file", "r");
}

void config_restore(toml_table_t *restore)
{
	if (restore == NULL)
		return;

	char *fileMode = "";

	toml_datum_t mode = toml_string_in(restore, "mode");
	if (!mode.ok)
	{
		error("Cannot read restore.mode", "Should be load / print / none");
	}

	if (strcmp(mode.u.s, "none") == 0)
	{
		cgen_config.load = NOTHING;

		free(mode.u.s);
		return;
	}
	else if (strcmp(mode.u.s, "print") == 0)
	{
		cgen_config.load = PRINT;
		fileMode = "w";
	}
	else if (strcmp(mode.u.s, "load") == 0)
	{
		cgen_config.load = LOAD_AND_PRINT;
		fileMode = "rw";
	}
	else
	{
		error("Unknown value for restore.mode", mode.u.s);
	}
	free(mode.u.s);

	cgen_config.currentFile = cgen_toml_file_in(restore, "current", fileMode);
	cgen_config.restartFile = cgen_toml_file_in(restore, "restart", fileMode);
}

static inline toml_datum_t cgen_toml_required(toml_datum_t method(const toml_table_t *, const char *), toml_table_t *t,
											  const char *key)
{
	toml_datum_t datum = method(t, key);
	if (!datum.ok)
		error("Key not found", getFullKey(t, key));

	return datum;
}

void config_output(toml_table_t *output)
{
	if (output == NULL)
		error("output section missing in TOML", "");

	cgen_config.energyFile = cgen_toml_file_in(output, "energy", "w");
	cgen_config.outputFile = cgen_toml_file_in(output, "coordinates", "w");

	toml_datum_t plots = toml_string_in(output, "plots");
	cgen_config.printPlot = plots.ok;
	if (!plots.ok)
	{
		fprintf(stderr, "output.plots not set, plot printing disabled. set output file to enable\n");
	}
	else
	{
		cgen_config.plotsFile = fopen(plots.u.s, "w");
		if (cgen_config.plotsFile == NULL)
		{
			error("Failed to open plots file for writting", plots.u.s);
		}
		free(plots.u.s);
	}
}

void config_simulation(toml_table_t *simulation)
{
	if (simulation == NULL)
		error("simulation section missing in TOML", "");

	toml_datum_t n = cgen_toml_required(toml_int_in, simulation, "N");
	cgen_config.N = n.u.i;

	toml_datum_t t0 = cgen_toml_required(toml_double_in, simulation, "T0");
	cgen_config.T0 = t0.u.d;

	toml_datum_t dt0 = cgen_toml_required(toml_double_in, simulation, "dt0");
	cgen_config.dt0 = dt0.u.d;

	toml_datum_t frames = cgen_toml_required(toml_int_in, simulation, "frames");
	cgen_config.frames = frames.u.i;

	toml_datum_t steps = cgen_toml_required(toml_int_in, simulation, "steps");
	cgen_config.steps = steps.u.i;

	toml_datum_t linear = cgen_toml_required(toml_bool_in, simulation, "linear");
	cgen_config.isLinear = linear.u.b;

	toml_datum_t seed = toml_int_in(simulation, "seed");
	if (!seed.ok)
	{
		cgen_config.seed = time(NULL);
		fprintf(stderr, "Cannot read simulation.seed, using random seed %ld\n", cgen_config.seed);
	}
	else
	{
		cgen_config.seed = seed.u.i;
	}
}

void config_dump(FILE *fd)
{
	fprintf(fd, "restore.mode = %d\n", cgen_config.load);

	// fprintf(fd, "restore.file = %s\n", cgen_config.restartFile);
	// fprintf(fd, "restore.current = %s\n", cgen_config.currentFile);
	//
	// fprintf(fd, "output.energy = %s\n", cgen_config.energyFile->_fileno);
	// fprintf(fd, "output.coordinates = %s\n", cgen_config.steps);

	fprintf(fd, "simulation.N = %d\n", cgen_config.N);
	fprintf(fd, "simulation.T0 = %lf\n", cgen_config.T0);
	fprintf(fd, "simulation.dt0 = %lf\n", cgen_config.dt0);
	fprintf(fd, "simulation.frames = %d\n", cgen_config.frames);
	fprintf(fd, "simulation.steps = %d\n", cgen_config.steps);
	fprintf(fd, "simulation.linear = %d\n", cgen_config.isLinear);
	fprintf(fd, "simulation.seed = %ld\n", cgen_config.seed);
}

#define ENV_PREFIX "CGEN_"
#define cgen_env_override_int(ENVNAME, VAR)                                                                            \
	do                                                                                                                 \
	{                                                                                                                  \
		char *var = getenv(ENV_PREFIX ENVNAME);                                                                        \
		if (var != NULL)                                                                                               \
		{                                                                                                              \
			errno = 0;                                                                                                 \
			VAR = strtoimax(var, NULL, 10);                                                                            \
			if (errno != 0)                                                                                            \
				error("Unsupported value in " ENV_PREFIX ENVNAME, var);                                                \
			fprintf(stderr, "WARN: using " ENVNAME "=%d\n", VAR);                                                      \
		}                                                                                                              \
	} while (0);

void config_process_env()
{
	cgen_env_override_int("FRAMES", cgen_config.frames);
	cgen_env_override_int("STEPS", cgen_config.steps);
}
