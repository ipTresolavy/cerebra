#include "bayes.h"

#include <math.h>

double crbr_log_sum_exp(const double *log_vals, int n)
{
        if(n <= 0)
        {
                return -INFINITY;
        }

        double max_val = log_vals[0];
        for(int i = 1; i < n; ++i)
        {
                if(log_vals[i] > max_val)
                {
                        max_val = log_vals[i];
                }
        }

        double sum = 0.0;
        for(int i = 0; i < n; ++i)
        {
                sum += exp(log_vals[i] - max_val);
        }

        return max_val + log(sum);
}

size_t crbr_bayes_classify(
        const double *log_priors, const double **log_likelihoods, size_t num_classes, size_t num_features
)
{
        size_t best_class    = -1;
        double best_log_prob = -INFINITY;

        for(size_t i = 0; i < num_classes; i++)
        {
                double log_prob = log_priors[i];
                for(size_t j = 0; j < num_features; j++)
                {
                        log_prob += log_likelihoods[i][j];
                }
                if(log_prob > best_log_prob)
                {
                        best_log_prob = log_prob;
                        best_class    = i;
                }
        }

        return best_class;
}
