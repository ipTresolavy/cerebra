#ifndef DENSITY_H
#define DENSITY_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

/*===========================================================================
 * Macros
 *===========================================================================*/

/**
 * @brief Macro definition for Pi.
 */
#ifndef CRBR_PI
        #define CRBR_PI M_PI
#endif

/*===========================================================================
 * Kernel Selection for KDE
 *===========================================================================*/

/**
 * @brief Enumeration for selecting the kernel function in KDE.
 */
typedef enum
{
        CRBR_KERNEL_GAUSSIAN,    /**< Gaussian kernel */
        CRBR_KERNEL_UNIFORM,     /**< Uniform kernel */
        CRBR_KERNEL_EPANECHNIKOV /**< Epanechnikov kernel */
} crbr_kernel_t;

/*===========================================================================
 * Gaussian Structures
 *===========================================================================*/

/**
 * @brief Structure representing a 1-dimensional Gaussian distribution.
 */
typedef struct crbr_gaussian_t
{
        double mean;     /**< Mean of the distribution */
        double variance; /**< Variance of the distribution */
} crbr_gaussian_t;

/**
 * @brief Structure representing a multivariate Gaussian distribution.
 *
 * The mean vector is of length `dim` and the covariance matrix is stored in row-major order,
 * with a total of dim x dim elements.
 *
 * This structure now caches the inverse of the covariance matrix and its determinant.
 */
typedef struct crbr_multivariate_gaussian_t
{
        size_t  dim;         /**< Dimensionality of the distribution */
        double *mean;        /**< Pointer to the mean vector (length dim) */
        double *cov;         /**< Pointer to the covariance matrix (dim x dim, row-major order) */
        bool    cov_inv_set; /**< Flag indicating if the inverse has been computed */
        double *cov_inv;     /**< Pointer to the cached inverse of the covariance matrix */
        double  cov_det;     /**< Determinant of the covariance matrix */
} crbr_multivariate_gaussian_t;

/*===========================================================================
 * 1-Dimensional Gaussian Density Functions
 *===========================================================================*/

/**
 * @brief Computes the Gaussian probability density function (PDF) for a given value.
 *
 * Uses the parameters from a 1-dimensional Gaussian structure.
 *
 * @param x The point at which to evaluate the density.
 * @param gaussian Pointer to a crbr_gaussian_t structure (contains mean and variance).
 * @return double The computed density at point x.
 */
double crbr_gaussian_pdf(double x, const crbr_gaussian_t *gaussian);

/**
 * @brief Computes the natural logarithm of the Gaussian probability density function.
 *
 * @param x The point at which to evaluate the log-density.
 * @param gaussian Pointer to a crbr_gaussian_t structure (contains mean and variance).
 * @return double The computed log-density at point x.
 */
double crbr_gaussian_log_pdf(double x, const crbr_gaussian_t *gaussian);

/**
 * @brief Estimates the probability density at a given point using Kernel Density Estimation (KDE).
 *
 * This version allows selection of the kernel type.
 *
 * @param data Pointer to an array of data points.
 * @param n Number of data points.
 * @param x The point at which to estimate the density.
 * @param bandwidth The bandwidth parameter.
 * @param kernel The kernel type to be used.
 * @return double The estimated density at point x.
 */
double crbr_kde(const double *data, size_t n, double x, double bandwidth, crbr_kernel_t kernel);

/*===========================================================================
 * Multidimensional Gaussian Density Functions
 *===========================================================================*/

/**
 * @brief Computes the multivariate Gaussian probability density function (PDF).
 *
 * Given a d-dimensional point x and the parameters of a multivariate Gaussian
 * (provided in the crbr_multivariate_gaussian_t structure), this function computes
 * the probability density.
 *
 * @param x Pointer to an array representing the d-dimensional point.
 * @param gaussian Pointer to a crbr_multivariate_gaussian_t structure.
 * @return double The computed density at point x.
 */
double crbr_multivariate_gaussian_pdf(const double *x, crbr_multivariate_gaussian_t *gaussian);

/**
 * @brief Computes the natural logarithm of the multivariate Gaussian probability density function.
 *
 * @param x Pointer to an array representing the d-dimensional point.
 * @param gaussian Pointer to a crbr_multivariate_gaussian_t structure.
 * @return double The computed log-density at point x.
 */
double crbr_multivariate_gaussian_log_pdf(const double *x, crbr_multivariate_gaussian_t *gaussian);

/*===========================================================================
 * Covariance Matrix Calculation
 *===========================================================================*/

/**
 * @brief Computes the sample covariance matrix from a set of data points.
 *
 * Given a dataset containing n samples of d-dimensional points (stored in row-major order, one sample per row),
 * this function computes the d x d sample covariance matrix.
 *
 * @param data Pointer to an array containing the data points (n x d).
 * @param n The number of samples.
 * @param dim The dimensionality (d) of each sample.
 * @param gaussian Pointer to a crbr_multivariate_gaussian_t structure where the covariance matrix will be stored (d x
 * d, row-major order).
 * @return int 0 on success, non-zero on error.
 */
int crbr_covariance_matrix(const double *data, size_t n, size_t dim, crbr_multivariate_gaussian_t *gaussian);

/*===========================================================================
 * Gaussian Maximum Likelihood Estimation
 *===========================================================================*/

/**
 * @brief Estimates the parameters of a 1-dimensional Gaussian via maximum likelihood.
 *
 * Given a set of data points, this function computes the maximum likelihood estimates
 * of the mean and variance, and stores them in the provided structure.
 *
 * @param data Pointer to an array of data points.
 * @param n The number of data points.
 * @param gaussian Pointer to a crbr_gaussian_t structure where the estimated parameters will be stored.
 * @return int 0 on success, non-zero on error.
 */
int crbr_fit_gaussian(const double *data, size_t n, crbr_gaussian_t *gaussian);

/**
 * @brief Estimates the parameters of a multivariate Gaussian via maximum likelihood.
 *
 * Given a dataset with n samples of d-dimensional points (stored in row-major order, one sample per row),
 * this function computes the maximum likelihood estimates of the mean vector and the covariance matrix.
 * Memory for the mean vector and covariance matrix will be allocated within the function and should be
 * freed by the user via crbr_free_multivariate_gaussian.
 *
 * Additionally, this function resets any previously cached inverse.
 *
 * @param data Pointer to an array containing the data points (n x d).
 * @param n The number of samples.
 * @param dim The dimensionality (d) of each sample.
 * @param gaussian Pointer to a crbr_multivariate_gaussian_t structure where the estimated parameters will be stored.
 * @return int 0 on success, non-zero on error.
 */
int crbr_fit_multivariate_gaussian(const double *data, size_t n, size_t dim, crbr_multivariate_gaussian_t *gaussian);

/**
 * @brief Frees the memory allocated for a multivariate Gaussian structure.
 *
 * This function frees the mean vector, covariance matrix, and cached covariance inverse stored in the given structure.
 *
 * @param gaussian Pointer to a crbr_multivariate_gaussian_t structure to be freed.
 */
void crbr_free_multivariate_gaussian(crbr_multivariate_gaussian_t *gaussian);

/*===========================================================================
 * Gaussian Mixture Model (GMM) Structures and Functions
 *===========================================================================*/

/**
 * @brief Structure representing a Gaussian Mixture Model (GMM).
 *
 * This structure contains the number of components, the dimensionality,
 * an array of component weights, and an array of multivariate Gaussian components.
 */
typedef struct
{
        size_t                        n_components; /**< Number of Gaussian components */
        size_t                        dim;          /**< Dimensionality of the data */
        double                       *weights;      /**< Array of component weights (length n_components) */
        crbr_multivariate_gaussian_t *components;   /**< Array of Gaussian components (length n_components) */
} crbr_gmm_t;

/**
 * @brief Fits a Gaussian Mixture Model (GMM) to the given data using the EM algorithm.
 *
 * @param data Pointer to an array of data points (n x d, row-major order).
 * @param n Number of data points.
 * @param dim Dimensionality (d) of each data point.
 * @param n_components Number of Gaussian components to fit.
 * @param max_iter Maximum number of EM iterations.
 * @param tol Convergence tolerance for the EM algorithm.
 * @param gmm Pointer to a crbr_gmm_t structure where the fitted model will be stored.
 * @return int 0 on success, non-zero on error.
 */
int crbr_fit_gmm(
        const double *data, size_t n, size_t dim, size_t n_components, size_t max_iter, double tol, crbr_gmm_t *gmm
);

/**
 * @brief Frees the memory allocated for a Gaussian Mixture Model.
 *
 * This function frees the memory associated with the component weights and Gaussian components.
 *
 * @param gmm Pointer to a crbr_gmm_t structure to be freed.
 */
void crbr_free_gmm(crbr_gmm_t *gmm);

#ifdef __cplusplus
}
#endif

#endif // DENSITY_H
