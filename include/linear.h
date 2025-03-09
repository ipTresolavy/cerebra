#ifndef LINEAR_H
#define LINEAR_H

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @brief Enumeration for specifying matrix memory layout.
 */
typedef enum
{
        CRBR_ROW_MAJOR, /**< Matrices are stored in row-major order */
        CRBR_COL_MAJOR  /**< Matrices are stored in column-major order */
} crbr_order_t;

/**
 * @brief Computes the dot product of two double vectors.
 *
 * @param n Number of elements in the vectors.
 * @param x Pointer to the first vector.
 * @param y Pointer to the second vector.
 * @return double The computed dot product.
 */
double crbr_dot(int n, const double *x, const double *y);

/**
 * @brief Performs matrix multiplication: C = A * B.
 *
 * @param order The memory layout of the matrices (CRBR_ROW_MAJOR or CRBR_COL_MAJOR).
 * @param m Number of rows in matrix A and C.
 * @param n Number of columns in matrix B and C.
 * @param k Number of columns in matrix A and rows in matrix B.
 * @param A Pointer to the first matrix.
 * @param B Pointer to the second matrix.
 * @param C Pointer to the output matrix.
 */
void crbr_gemm(crbr_order_t order, int m, int n, int k, const double *A, const double *B, double *C);

/**
 * @brief Scales a vector by a constant factor: x = alpha * x.
 *
 * @param n Number of elements in the vector.
 * @param alpha Scaling factor.
 * @param x Pointer to the vector to scale.
 * @param incx Stride within vector x.
 */
void crbr_scal(int n, double alpha, double *x, int incx);

/**
 * @brief Computes the axpy operation: y = a*x + y.
 *
 * @param n Number of elements in the vectors.
 * @param a Scalar multiplier.
 * @param x Pointer to the input vector x.
 * @param y Pointer to the input/output vector y.
 * @param incx Stride within vector x.
 * @param incy Stride within vector y.
 */
void crbr_axpy(int n, double a, const double *x, double *y, int incx, int incy);

/**
 * @brief Performs matrix-vector multiplication: y = A * x.
 *
 * @param order The memory layout of the matrix A.
 * @param m Number of rows in matrix A.
 * @param n Number of columns in matrix A.
 * @param A Pointer to the matrix A.
 * @param x Pointer to the vector x.
 * @param y Pointer to the result vector y.
 */
void crbr_gemv(crbr_order_t order, int m, int n, const double *A, const double *x, double *y);

/**
 * @brief Performs the rank-1 update: A = A + alpha * x * y^T.
 *
 * @param order The memory layout of the matrix A.
 * @param m Number of rows in matrix A.
 * @param n Number of columns in matrix A.
 * @param alpha Scalar multiplier.
 * @param x Pointer to the vector x.
 * @param incx Stride within vector x.
 * @param y Pointer to the vector y.
 * @param incy Stride within vector y.
 * @param A Pointer to the matrix A to be updated.
 */
void crbr_ger(
        crbr_order_t order, int m, int n, double alpha, const double *x, int incx, const double *y, int incy, double *A
);

#ifdef __cplusplus
}
#endif

#endif // LINEAR_H
