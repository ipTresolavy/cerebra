#include "density.h"

#include "linear.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*==============================================================================
 * Helper Function: Gauss–Jordan Inversion with Determinant Calculation
 *============================================================================*/

/**
 * @brief Computes the inverse and determinant of a d x d matrix using Gauss–Jordan elimination.
 *
 * This helper is used to compute the inverse and determinant of the covariance matrix.
 */
static int gauss_jordan_inverse_det(const double *A, double *A_inv, size_t d, double *det)
{
        size_t  i, j, k;
        double *M = malloc(d * d * sizeof(double));
        if(M == NULL)
        {
                return -1;
        }
        memcpy(M, A, d * d * sizeof(double));

        // Initialize A_inv as the identity matrix.
        for(i = 0; i < d; i++)
        {
                for(j = 0; j < d; j++)
                {
                        A_inv[i * d + j] = (i == j) ? 1.0 : 0.0;
                }
        }

        *det = 1.0;
        for(i = 0; i < d; i++)
        {
                size_t pivot   = i;
                double max_val = fabs(M[i * d + i]);
                for(j = i + 1; j < d; j++)
                {
                        double val = fabs(M[j * d + i]);
                        if(val > max_val)
                        {
                                max_val = val;
                                pivot   = j;
                        }
                }
                if(fabs(M[pivot * d + i]) < 1e-12)
                { // Singular matrix
                        free(M);
                        return -1;
                }
                if(pivot != i)
                {
                        for(j = 0; j < d; j++)
                        {
                                double tmp       = M[i * d + j];
                                M[i * d + j]     = M[pivot * d + j];
                                M[pivot * d + j] = tmp;

                                tmp                  = A_inv[i * d + j];
                                A_inv[i * d + j]     = A_inv[pivot * d + j];
                                A_inv[pivot * d + j] = tmp;
                        }
                        *det = -(*det);
                }
                double pivot_val  = M[i * d + i];
                *det             *= pivot_val;
                for(j = 0; j < d; j++)
                {
                        M[i * d + j]     /= pivot_val;
                        A_inv[i * d + j] /= pivot_val;
                }
                for(j = 0; j < d; j++)
                {
                        if(j != i)
                        {
                                double factor = M[j * d + i];
                                for(k = 0; k < d; k++)
                                {
                                        M[j * d + k]     -= factor * M[i * d + k];
                                        A_inv[j * d + k] -= factor * A_inv[i * d + k];
                                }
                        }
                }
        }
        free(M);
        return 0;
}

/*==============================================================================
 * 1-Dimensional Gaussian Functions
 *============================================================================*/

double crbr_gaussian_pdf(double x, const crbr_gaussian_t *gaussian)
{
        if(NULL == gaussian)
        {
                return 0.0;
        }
        double mean     = gaussian->mean;
        double variance = gaussian->variance;
        double denom    = sqrt(2 * CRBR_PI * variance);
        double exponent = -((x - mean) * (x - mean)) / (2 * variance);
        return (denom < 1e-12) ? 0.0 : (exp(exponent) / denom);
}

double crbr_gaussian_log_pdf(double x, const crbr_gaussian_t *gaussian)
{
        if(NULL == gaussian)
        {
                return 0.0;
        }
        double mean     = gaussian->mean;
        double variance = gaussian->variance;
        if(variance < 1e-12)
        {
                return 0.0;
        }
        return -0.5 * log(2 * CRBR_PI * variance) - ((x - mean) * (x - mean)) / (2 * variance);
}

/*==============================================================================
 * Multidimensional Gaussian Functions
 *============================================================================*/

double crbr_multivariate_gaussian_pdf(const double *x, crbr_multivariate_gaussian_t *gaussian)
{
        if(NULL == x || NULL == gaussian)
        {
                return 0.0;
        }
        size_t        dim  = gaussian->dim;
        const double *mean = gaussian->mean;
        const double *cov  = gaussian->cov;

        double  det;
        double *cov_inv = NULL;
        if(gaussian->cov_inv_set)
        {
                cov_inv = gaussian->cov_inv;
                det     = gaussian->cov_det;
        }
        else
        {
                cov_inv = malloc(dim * dim * sizeof(double));
                if(cov_inv == NULL)
                {
                        return 0.0;
                }
                if(gauss_jordan_inverse_det(cov, cov_inv, dim, &det) != 0)
                {
                        free(cov_inv);
                        return 0.0;
                }
                gaussian->cov_inv     = cov_inv;
                gaussian->cov_inv_set = true;
                gaussian->cov_det     = det;
        }

        double *diff = malloc(dim * sizeof(double));
        if(diff == NULL)
        {
                return 0.0;
        }
        for(size_t i = 0; i < dim; i++)
        {
                diff[i] = x[i] - mean[i];
        }

        double *temp = malloc(dim * sizeof(double));
        if(temp == NULL)
        {
                free(diff);
                return 0.0;
        }
        crbr_gemv(CRBR_ROW_MAJOR, (int)dim, (int)dim, cov_inv, diff, temp);

        double mahal = crbr_dot((int)dim, diff, temp);

        free(diff);
        free(temp);

        double norm_const = pow(2 * CRBR_PI, (double)dim / 2.0) * sqrt(det);
        return exp(-0.5 * mahal) / norm_const;
}

double crbr_multivariate_gaussian_log_pdf(const double *x, crbr_multivariate_gaussian_t *gaussian)
{
        if(NULL == x || NULL == gaussian)
        {
                return -INFINITY;
        }
        size_t        dim  = gaussian->dim;
        const double *mean = gaussian->mean;
        const double *cov  = gaussian->cov;

        double  det;
        double *cov_inv = NULL;
        if(gaussian->cov_inv_set)
        {
                cov_inv = gaussian->cov_inv;
                det     = gaussian->cov_det;
        }
        else
        {
                cov_inv = malloc(dim * dim * sizeof(double));
                if(cov_inv == NULL)
                {
                        return -INFINITY;
                }
                if(gauss_jordan_inverse_det(cov, cov_inv, dim, &det) != 0)
                {
                        free(cov_inv);
                        return -INFINITY;
                }
                gaussian->cov_inv     = cov_inv;
                gaussian->cov_inv_set = true;
                gaussian->cov_det     = det;
        }

        double *diff = malloc(dim * sizeof(double));
        if(diff == NULL)
        {
                return -INFINITY;
        }
        for(size_t i = 0; i < dim; i++)
        {
                diff[i] = x[i] - mean[i];
        }

        double *temp = malloc(dim * sizeof(double));
        if(temp == NULL)
        {
                free(diff);
                return -INFINITY;
        }
        crbr_gemv(CRBR_ROW_MAJOR, (int)dim, (int)dim, cov_inv, diff, temp);

        double mahal = crbr_dot((int)dim, diff, temp);

        free(diff);
        free(temp);

        double log_norm = ((double)dim / 2.0) * log(2 * CRBR_PI) + 0.5 * log(det);
        return -0.5 * mahal - log_norm;
}

/*==============================================================================
 * Covariance Matrix Calculation (Updated for Cached Inverse)
 *============================================================================*/

int crbr_covariance_matrix(const double *data, size_t n, size_t dim, crbr_multivariate_gaussian_t *gaussian)
{
        if(n < 2 || NULL == data || NULL == gaussian)
        {
                return -1;
        }
        // Allocate the mean vector.
        gaussian->mean = calloc(dim, sizeof(double));
        if(gaussian->mean == NULL)
        {
                return -1;
        }
        // Allocate the covariance matrix.
        gaussian->cov = calloc(dim * dim, sizeof(double));
        if(gaussian->cov == NULL)
        {
                free(gaussian->mean);
                gaussian->mean = NULL;
                return -1;
        }
        // Compute mean vector.
        for(size_t i = 0; i < n; i++)
        {
                for(size_t j = 0; j < dim; j++)
                {
                        gaussian->mean[j] += data[i * dim + j];
                }
        }
        for(size_t j = 0; j < dim; j++)
        {
                gaussian->mean[j] /= n;
        }
        // Zero out the covariance matrix.
        memset(gaussian->cov, 0, dim * dim * sizeof(double));
        double *diff = malloc(dim * sizeof(double));
        if(diff == NULL)
        {
                free(gaussian->mean);
                free(gaussian->cov);
                gaussian->mean = NULL;
                gaussian->cov  = NULL;
                return -1;
        }
        // For each sample, compute diff = (x - mean) and update covariance via a rank-1 update.
        for(size_t i = 0; i < n; i++)
        {
                for(size_t j = 0; j < dim; j++)
                {
                        diff[j] = data[i * dim + j] - gaussian->mean[j];
                }
                crbr_ger(CRBR_ROW_MAJOR, (int)dim, (int)dim, 1.0, diff, 1, diff, 1, gaussian->cov);
        }
        // Divide the accumulated sum by (n - 1)
        for(size_t i = 0; i < dim * dim; i++)
        {
                gaussian->cov[i] /= (n - 1);
        }
        free(diff);
        // Reset any previously cached inverse.
        if(gaussian->cov_inv_set)
        {
                if(gaussian->cov_inv != NULL)
                {
                        free(gaussian->cov_inv);
                }
                gaussian->cov_inv = NULL;
        }
        gaussian->cov_inv_set = false;
        gaussian->cov_det     = 0.0;
        return 0;
}

/*==============================================================================
 * Maximum Likelihood Estimation for Gaussian Parameters
 *============================================================================*/

int crbr_fit_gaussian(const double *data, size_t n, crbr_gaussian_t *gaussian)
{
        if(n == 0 || NULL == data || NULL == gaussian)
        {
                return -1;
        }
        double sum = 0.0;
        for(size_t i = 0; i < n; i++)
        {
                sum += data[i];
        }
        double mean    = sum / n;
        double var_sum = 0.0;
        for(size_t i = 0; i < n; i++)
        {
                double diff  = data[i] - mean;
                var_sum     += diff * diff;
        }
        gaussian->mean     = mean;
        gaussian->variance = var_sum / n; // ML estimator uses 1/n
        return 0;
}

int crbr_fit_multivariate_gaussian(const double *data, size_t n, size_t dim, crbr_multivariate_gaussian_t *gaussian)
{
        if(n == 0 || NULL == data || NULL == gaussian)
        {
                return -1;
        }
        gaussian->dim  = dim;
        gaussian->mean = calloc(dim, sizeof(double));
        gaussian->cov  = calloc(dim * dim, sizeof(double));
        if(gaussian->mean == NULL || gaussian->cov == NULL)
        {
                free(gaussian->mean);
                free(gaussian->cov);
                return -1;
        }
        for(size_t i = 0; i < n; i++)
        {
                for(size_t j = 0; j < dim; j++)
                {
                        gaussian->mean[j] += data[i * dim + j];
                }
        }
        for(size_t j = 0; j < dim; j++)
        {
                gaussian->mean[j] /= n;
        }
        memset(gaussian->cov, 0, dim * dim * sizeof(double));
        double *diff = malloc(dim * sizeof(double));
        if(diff == NULL)
        {
                return -1;
        }
        for(size_t i = 0; i < n; i++)
        {
                for(size_t j = 0; j < dim; j++)
                {
                        diff[j] = data[i * dim + j] - gaussian->mean[j];
                }
                crbr_ger(CRBR_ROW_MAJOR, (int)dim, (int)dim, 1.0, diff, 1, diff, 1, gaussian->cov);
        }
        for(size_t i = 0; i < dim * dim; i++)
        {
                gaussian->cov[i] /= n;
        }
        free(diff);
        // Reset any previously cached inverse.
        gaussian->cov_inv_set = false;
        if(gaussian->cov_inv != NULL)
        {
                free(gaussian->cov_inv);
                gaussian->cov_inv = NULL;
        }
        gaussian->cov_det = 0.0;
        return 0;
}

void crbr_free_multivariate_gaussian(crbr_multivariate_gaussian_t *gaussian)
{
        if(gaussian)
        {
                free(gaussian->mean);
                free(gaussian->cov);
                if(gaussian->cov_inv != NULL)
                {
                        free(gaussian->cov_inv);
                }
                gaussian->mean        = NULL;
                gaussian->cov         = NULL;
                gaussian->cov_inv     = NULL;
                gaussian->cov_inv_set = false;
        }
}

/*==============================================================================
 * Gaussian Mixture Model (GMM) via EM Algorithm
 *============================================================================*/

int crbr_fit_gmm(
        const double *data, size_t n, size_t dim, size_t n_components, size_t max_iter, double tol, crbr_gmm_t *gmm
)
{
        if(n == 0 || n_components == 0 || NULL == data || NULL == gmm)
        {
                return -1;
        }
        gmm->n_components = n_components;
        gmm->dim          = dim;
        gmm->weights      = calloc(n_components, sizeof(double));
        gmm->components   = calloc(n_components, sizeof(crbr_multivariate_gaussian_t));
        if(gmm->weights == NULL || gmm->components == NULL)
        {
                free(gmm->weights);
                free(gmm->components);
                return -1;
        }
        // Initialize weights uniformly.
        for(size_t i = 0; i < n_components; i++)
        {
                gmm->weights[i] = 1.0 / n_components;
        }
        // Initialize each component with a data point and identity covariance.
        for(size_t i = 0; i < n_components; i++)
        {
                gmm->components[i].dim  = dim;
                gmm->components[i].mean = calloc(dim, sizeof(double));
                gmm->components[i].cov  = calloc(dim * dim, sizeof(double));
                if(gmm->components[i].mean == NULL || gmm->components[i].cov == NULL)
                {
                        for(size_t j = 0; j < i; j++)
                        {
                                crbr_free_multivariate_gaussian(&gmm->components[j]);
                        }
                        free(gmm->weights);
                        free(gmm->components);
                        return -1;
                }
                size_t idx = (i < n) ? i : (n - 1);
                for(size_t j = 0; j < dim; j++)
                {
                        gmm->components[i].mean[j] = data[idx * dim + j];
                }
                for(size_t j = 0; j < dim; j++)
                {
                        for(size_t k = 0; k < dim; k++)
                        {
                                gmm->components[i].cov[j * dim + k] = (j == k) ? 1.0 : 0.0;
                        }
                }
        }
        double *resp = calloc(n * n_components, sizeof(double));
        if(resp == NULL)
        {
                for(size_t i = 0; i < n_components; i++)
                {
                        crbr_free_multivariate_gaussian(&gmm->components[i]);
                }
                free(gmm->weights);
                free(gmm->components);
                return -1;
        }
        double prev_log_likelihood = -INFINITY;
        for(size_t iter = 0; iter < max_iter; iter++)
        {
                double log_likelihood = 0.0;
                // E-step: compute responsibilities.
                for(size_t i = 0; i < n; i++)
                {
                        double sum_resp = 0.0;
                        for(size_t j = 0; j < n_components; j++)
                        {
                                double pdf = crbr_multivariate_gaussian_pdf(&data[i * dim], &gmm->components[j]);
                                double r   = gmm->weights[j] * pdf;
                                resp[i * n_components + j]  = r;
                                sum_resp                   += r;
                        }
                        if(sum_resp == 0)
                        {
                                sum_resp = 1e-12;
                        }
                        log_likelihood += log(sum_resp);
                        for(size_t j = 0; j < n_components; j++)
                        {
                                resp[i * n_components + j] /= sum_resp;
                        }
                }
                if(fabs(log_likelihood - prev_log_likelihood) < tol)
                {
                        break;
                }
                prev_log_likelihood = log_likelihood;
                // M-step: update parameters.
                for(size_t j = 0; j < n_components; j++)
                {
                        double resp_sum = 0.0;
                        for(size_t i = 0; i < n; i++)
                        {
                                resp_sum += resp[i * n_components + j];
                        }
                        gmm->weights[j] = resp_sum / n;
                        // Update mean.
                        for(size_t k = 0; k < dim; k++)
                        {
                                double weighted_sum = 0.0;
                                for(size_t i = 0; i < n; i++)
                                {
                                        weighted_sum += resp[i * n_components + j] * data[i * dim + k];
                                }
                                if(resp_sum > 0)
                                {
                                        gmm->components[j].mean[k] = weighted_sum / resp_sum;
                                }
                        }
                        // Update covariance.
                        memset(gmm->components[j].cov, 0, dim * dim * sizeof(double));
                        double *diff = malloc(dim * sizeof(double));
                        if(diff == NULL)
                        {
                                continue;
                        }
                        for(size_t i = 0; i < n; i++)
                        {
                                for(size_t k = 0; k < dim; k++)
                                {
                                        diff[k] = data[i * dim + k] - gmm->components[j].mean[k];
                                }
                                crbr_ger(
                                        CRBR_ROW_MAJOR, (int)dim, (int)dim, resp[i * n_components + j], diff, 1, diff,
                                        1, gmm->components[j].cov
                                );
                        }
                        free(diff);
                        for(size_t k = 0; k < dim * dim; k++)
                        {
                                gmm->components[j].cov[k] /= resp_sum;
                        }
                }
        }
        free(resp);
        return 0;
}

void crbr_free_gmm(crbr_gmm_t *gmm)
{
        if(gmm)
        {
                if(gmm->weights)
                {
                        free(gmm->weights);
                        gmm->weights = NULL;
                }
                if(gmm->components)
                {
                        for(size_t i = 0; i < gmm->n_components; i++)
                        {
                                crbr_free_multivariate_gaussian(&gmm->components[i]);
                        }
                        free(gmm->components);
                        gmm->components = NULL;
                }
        }
}
