#include <stdio.h>
#include <stdlib.h>
#define COMPEARTH_PRIVATE_DET3X3 1
#define COMPEARTH_PRIVATE_GEMV3 1
#define COMPEARTH_PRIVATE_GEM3 1
#define COMPEARTH_PRIVATE_GEMT3 1
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_NORM3 1
#define COMPEARTH_PRIVATE_DOT3 1
#define COMPEARTH_PRIVATE_WRAP360 1
#include "compearth.h"

double det3x3ColumnMajor(const double *__restrict__ A)
{
    double det;
    det = A[0]*( A[4]*A[8] - A[5]*A[7]) 
        - A[3]*( A[1]*A[8] - A[2]*A[7])
        + A[6]*( A[1]*A[5] - A[2]*A[4]);
    return det;
};

void gemm3_colMajorNoTransNoTrans(const double *__restrict__ A,
                                  const double *__restrict__ B,
                                  double *__restrict__ C)
{
    // column 1
    C[0] = A[0]*B[0] + A[3]*B[1] + A[6]*B[2];
    C[1] = A[1]*B[0] + A[4]*B[1] + A[7]*B[2];
    C[2] = A[2]*B[0] + A[5]*B[1] + A[8]*B[2];
    // column 2
    C[3] = A[0]*B[3] + A[3]*B[4] + A[6]*B[5];
    C[4] = A[1]*B[3] + A[4]*B[4] + A[7]*B[5];
    C[5] = A[2]*B[3] + A[5]*B[4] + A[8]*B[5];
    // column 3
    C[6] = A[0]*B[6] + A[3]*B[7] + A[6]*B[8];
    C[7] = A[1]*B[6] + A[4]*B[7] + A[7]*B[8];
    C[8] = A[2]*B[6] + A[5]*B[7] + A[8]*B[8];
    return; 
}

void gemm3_colMajorNoTransTrans(const double *__restrict__ A,
                                const double *__restrict__ B,
                                double *__restrict__ C)
{
    // column 1
    C[0] = A[0]*B[0] + A[3]*B[3] + A[6]*B[6];
    C[1] = A[1]*B[0] + A[4]*B[3] + A[7]*B[6];
    C[2] = A[2]*B[0] + A[5]*B[3] + A[8]*B[6];
    // column 2
    C[3] = A[0]*B[1] + A[3]*B[4] + A[6]*B[7];
    C[4] = A[1]*B[1] + A[4]*B[4] + A[7]*B[7];
    C[5] = A[2]*B[1] + A[5]*B[4] + A[8]*B[7];
    // column 3
    C[6] = A[0]*B[2] + A[3]*B[5] + A[6]*B[8];
    C[7] = A[1]*B[2] + A[4]*B[5] + A[7]*B[8];
    C[8] = A[2]*B[2] + A[5]*B[5] + A[8]*B[8];
    return;  
}

void gemv3_colMajorNoTrans(const double *__restrict__ A,
                           const double *__restrict__ x,
                           double *__restrict__ y)
{
    y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
    y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
    y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
    return;
}

void cross3(const double *__restrict__ a,  
            const double *__restrict__ b,
            double *__restrict__ c)
{
    // Compute cross product
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

double norm3(const double *__restrict__ a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

double dot3(const double *__restrict__ a, const double *__restrict__ b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/*!
 * @brief Wraps angle in degrees to [0,360]
 *
 * @param[in] lon   angle to wrap (degrees)
 *
 * @result wrapped angle in range [0,360]
 *
 * @author Carl Tape translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
double wrap360(const double lon)
{
    double lonw;
    bool lpos;
    lpos = false; 
    if (lon > 0.0){lpos = true;}
    lonw = fmod(lon, 360.0); 
    if (lonw == 0.0 && lpos){lonw = 360.0;}
    return lonw;
}
