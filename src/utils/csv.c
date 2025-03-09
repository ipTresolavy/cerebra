#include "utils/csv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *load_csv(const char *filename, size_t *out_rows, size_t *out_cols)
{
        FILE *fp = fopen(filename, "r");
        if(fp == NULL)
        {
                fprintf(stderr, "Error opening file %s\n", filename);
                return NULL;
        }

        char   line[4096];
        size_t rows = 0;
        size_t cols = 0;

        // First pass: determine number of rows and columns.
        if(fgets(line, sizeof(line), fp) != NULL)
        {
                rows++;
                char *token = strtok(line, ",");
                while(token != NULL)
                {
                        cols++;
                        token = strtok(NULL, ",");
                }
        }
        while(fgets(line, sizeof(line), fp) != NULL)
        {
                rows++;
        }

        // Allocate the data array.
        double *data = malloc(rows * cols * sizeof(double));
        if(data == NULL)
        {
                fclose(fp);
                return NULL;
        }

        // Rewind file for second pass.
        fseek(fp, 0, SEEK_SET);
        size_t row = 0;
        while(fgets(line, sizeof(line), fp) != NULL)
        {
                size_t col   = 0;
                char  *token = strtok(line, ",");
                while(token != NULL && col < cols)
                {
                        data[row * cols + col] = strtod(token, NULL);
                        col++;
                        token = strtok(NULL, ",");
                }
                row++;
        }

        fclose(fp);
        *out_rows = rows;
        *out_cols = cols;
        return data;
}
