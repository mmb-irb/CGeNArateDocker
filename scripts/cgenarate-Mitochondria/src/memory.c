
#include <stdio.h>
#include <stdlib.h>

FILE *fopen_read(char *name, char *error_msg)
{
	FILE *file;
	file = fopen(name, "r");
	if (file == NULL)
	{
		printf("%s", error_msg);
		exit(1);
	}

	return file;
}

FILE *fopen_write(char *name, char *error_msg)
{
	FILE *file;
	file = fopen(name, "w");
	if (file == NULL)
	{
		printf("%s", error_msg);
		exit(1);
	}

	return file;
}

void PrintCoordinates(double (*_r)[][3], int n, FILE *coordinatesfile)
{
	double(*r)[n][3] = _r;
	int i, l;
	for (i = 0; i < n; i++)
		for (l = 0; l < 3; l++)
		{
			fprintf(coordinatesfile, "%8.2lf", r[0][i][l]);
			if ((3 * i + l + 1) % 10 == 0)
				fprintf(coordinatesfile, "\n");
		}
	if ((3 * (i - 1) + l) % 10 != 0)
		fprintf(coordinatesfile, "\n");
}
