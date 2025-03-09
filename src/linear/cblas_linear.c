#include "linear.h"

#include <cblas/cblas.h>

double crbr_dot(int n, const double *x, const double *y)
{
        return cblas_ddot(n, x, 1, y, 1);
}

void crbr_gemm(crbr_order_t order, int m, int n, int k, const double *A, const double *B, double *C)
{
        if(order == CRBR_ROW_MAJOR)
        {
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);
        }
        else
        { // CRBR_COL_MAJOR
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, m, B, k, 0.0, C, m);
        }
}

void crbr_scal(int n, double alpha, double *x, int incx)
{
        cblas_dscal(n, alpha, x, incx);
}

void crbr_axpy(int n, double a, const double *x, double *y, int incx, int incy)
{
        cblas_daxpy(n, a, x, incx, y, incy);
}

void crbr_gemv(crbr_order_t order, int m, int n, const double *A, const double *x, double *y)
{
        if(order == CRBR_ROW_MAJOR)
        {
                cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, 1, 0.0, y, 1);
        }
        else
        {
                cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, A, m, x, 1, 0.0, y, 1);
        }
}

void crbr_ger(
        crbr_order_t order, int m, int n, double alpha, const double *x, int incx, const double *y, int incy, double *A
)
{
        if(order == CRBR_ROW_MAJOR)
        {
                cblas_dger(CblasRowMajor, m, n, alpha, x, incx, y, incy, A, n);
        }
        else
        {
                cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, m);
        }
}
