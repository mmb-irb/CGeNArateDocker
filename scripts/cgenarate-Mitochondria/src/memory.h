#ifndef MEMORY_H
#define MEMORY_H

#include <stdio.h>

FILE *fopen_read(char *name, char *error_msg);
FILE *fopen_write(char *name, char *error_msg);

void PrintCoordinates(double (*)[][3], int n, FILE *coordinatesfile);

#endif // MEMORY_H
