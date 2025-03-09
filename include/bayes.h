#ifndef BAYES_H
#define BAYES_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>

/**
 * @brief Computes the log-sum-exp of an array of log values.
 *
 * This function safely computes the logarithm of the sum of exponentials
 * of the input values, a common operation in Bayes decision theory to avoid
 * numerical underflow.
 *
 * @param log_vals Pointer to an array of log values.
 * @param n The number of elements in the array.
 * @return double The computed log-sum-exp value.
 */
double crbr_log_sum_exp(const double *log_vals, int n);

/**
 * @brief Classifies an instance using Bayes decision theory.
 *
 * This function calculates the log-posterior for each class by summing the
 * log prior and the log likelihoods for each feature, and then returns the
 * index of the class with the highest log probability.
 *
 * @param log_priors Pointer to an array containing the log prior probabilities for each class.
 * @param log_likelihoods A pointer to an array of pointers, where each sub-array contains the log likelihoods for the
 * features corresponding to a class.
 * @param num_classes The number of classes.
 * @param num_features The number of features.
 * @return int The index of the class with the highest posterior probability.
 */
size_t crbr_bayes_classify(
        const double *log_priors, const double **log_likelihoods, size_t num_classes, size_t num_features
);

#ifdef __cplusplus
}
#endif

#endif // BAYES_H
