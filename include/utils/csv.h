#ifndef CSV_H
#define CSV_H

#include <stddef.h>

/**
 * @brief Loads a CSV file with numerical data.
 *
 * The CSV file is assumed to have no header and to use commas as delimiters.
 *
 * @param filename The name of the CSV file to load.
 * @param out_rows Pointer to a size_t where the number of rows will be stored.
 * @param out_cols Pointer to a size_t where the number of columns will be stored.
 * @return double* Pointer to the allocated array of doubles in row-major order, or NULL on error.
 */
double *load_csv(const char *filename, size_t *out_rows, size_t *out_cols);

#endif // CSV_H
