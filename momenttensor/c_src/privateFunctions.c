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
#define COMPEARTH_PRIVATE_MOD 1
#define COMPEARTH_PRIVATE_ANTIPODE 1
#define COMPEARTH_PRIVATE_UPDOWN_ABS_ARGSORT3 1
#define COMPEARTH_PRIVATE_UPDOWN_ARGSORT3 1
#define COMPEARTH_PRIVATE_ARGSORT3 1
#include "compearth.h"


/*!
 * @brief Computes the determinant of a 3 x 3 major in column major order.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline double det3x3ColumnMajor(const double *__restrict__ A)
{
    double det;
    det = A[0]*( A[4]*A[8] - A[5]*A[7]) 
        - A[3]*( A[1]*A[8] - A[2]*A[7])
        + A[6]*( A[1]*A[5] - A[2]*A[4]);
    return det;
}
/*!
 * @brief Computes C = A*B where A, B, and C are 3 x 3 matrices in column
 *        major order.
 */ 
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline void gemm3_colMajorNoTransNoTrans(const double *__restrict__ A,
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
/*!
 * @brief Computes C = A*B' where A, B, and C are 3 x 3 matrices in column
 *        major order.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline void gemm3_colMajorNoTransTrans(const double *__restrict__ A,
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
/*!
 * @brief Computes y = A*x where A is a 3 x 3 matrix in column major order
 *        and x and y are length 3 vectors.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline void gemv3_colMajorNoTrans(const double *__restrict__ A,
                                  const double *__restrict__ x,
                                  double *__restrict__ y)
{
    y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
    y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
    y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
    return;
}
/*!
 * @brief Computes the cross-product, c = a x c, where a, b, and c are 
 *        length 3 vectors.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline void cross3(const double *__restrict__ a,  
                   const double *__restrict__ b,
                   double *__restrict__ c)
{
    // Compute cross product
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}
/*!
 * @brief Computes the norm of a vector a which is length 3.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline double norm3(const double *__restrict__ a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
/*!
 * @brief Computes the dot-product a.b where a and b are length 3 vectors.
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline double dot3(const double *__restrict__ a, const double *__restrict__ b)
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
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline double wrap360(const double lon)
{
    double lonw;
    bool lpos;
    lpos = false; 
    if (lon > 0.0){lpos = true;}
    lonw = lon - floor(lon/360.0)*360.0; //lonw = fmod(lon, 360.0); 
    if (lon == 360.){lonw = 360.0;} // matlab convention
    if (lonw == 0.0 && lpos){lonw = 360.0;}
    return lonw;
}

/*!
 * @brief Emulation of Matlab's modulus after division.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
#ifdef _OPENMP
#pragma omp declare simd
#endif
inline double mod(const double x, const double y)
{
    double xmod;
    bool xisnty;
    xmod = 0.0; // Convention 1 - if y == 0.0
    // This will avoid division by 0
    if (fabs(y) > 0.0)
    {
        xmod = x - floor(x/y)*y;
        // Convention 2 - mod(x, x) is 0
        xisnty = false;
        if (x != y){xisnty = true;}
        if (!xisnty){xmod = 0.0;} //if (x == y){xmod = 0.0;}
        // Convention 3 - if x ~= y and y~=0 then mod has the same
        // sign as y.  Note y~= 0 is true b/c fabs(y) > 0.
        if (xisnty && y < 0.0 && xmod > 0.0){xmod =-xmod;}
        if (xisnty && y > 0.0 && xmod < 0.0){xmod =+xmod;}
    }
/*
    // Original code
    if (y == 0.0){return x;}
    if (x == y){return 0.0;}
    xmod = x - floor(x/y)*y;
    // 3rd convention - for x ~= y and y ~= 0 mod has same sign as y
    if (x != y && y != 0.0)
    {
        if (y < 0.0){xmod =-xmod;}
        if (y > 0.0){xmod =+xmod;}
    }
*/
    return xmod; 
}

#ifdef _OPENMP
#pragma omp declare simd
#endif
inline void antipode(const double lat, const double lon,
                     double *latOut, double *lonOut, const bool isDeg)
{
    *latOut =-lat;
    if (isDeg)
    {
        *lonOut = 180.0 - mod(-lon, 360.0);
    }
    else
    {
        *lonOut = M_PI - mod(-lon, 2.0*M_PI);
    }
    return;
}

/*!
 * @brief Argsort a length 3 vector in increasing or decreasing order.
 *
 * @param[in] x         Length 3 array to argsort.
 * @param[in] lascend   If true then sort in ascending order. \n
 *                      Otherwise, sort in descending order.
 *
 * @param[out] perm     Permutation such that x(perm) is in ascending
 *                      or descending order.
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 * @copyright MIT
 *
 */
inline int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ perm)
{
    int permt[3], ierr;
    if (lascend)
    {
        ierr = argsort3(x, perm);
    }
    else
    {
        ierr = argsort3(x, permt);
        perm[0] = permt[2];
        perm[1] = permt[1];
        perm[2] = permt[0];
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Argsort a length 3 vector in increasing or decreasing order 
 *        of absolute value.
 *
 * @param[in] x         Length 3 array to argsort.
 * @param[in] lascend   If true then sort in ascending order. \n
 *                      Otherwise, sort in descending order.
 *
 * @param[out] perm     Permutation such that |x(perm)| is in ascending
 *                      or descending order.
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 * @copyright MIT
 *
 */
inline int argsort3_absUpDown(const double *__restrict__ x,
                              const bool lascend,
                              int *__restrict__ perm)
{
    double xa[3];
    int ierr;
    xa[0] = fabs(x[0]);
    xa[1] = fabs(x[1]);
    xa[2] = fabs(x[2]);
    ierr = argsort3_upDown(xa, lascend, perm);
    return ierr;
}
//============================================================================//
/*!
 * @brief Argsorts a length 3 array in ascending order.
 *
 * @param[in] x      Length 3 array to argsort.
 *
 * @param[out] perm  Permutation such that x(perm) is in ascending order.
 * 
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 * @copyright MIT
 *
 */ 
inline int argsort3(const double *__restrict__ x, int *__restrict__ perm)
{
    int i, temp;
    const int a = 0;
    const int b = 1;
    const int c = 2;
    // Copy
    perm[a] = a;
    perm[b] = b;
    perm[c] = c;
    if (x[perm[a]] > x[perm[c]])
    {
        temp = perm[c];
        perm[c] = perm[a];
        perm[a] = temp;
    }
    if (x[perm[a]] > x[perm[b]])
    {
        temp = perm[b];
        perm[b] = perm[a];
        perm[a] = temp;
    }
    //Now the smallest element is the first one. Just check the 2-nd and 3-rd
    if (x[perm[b]] > x[perm[c]])
    {
        temp = perm[c];
        perm[c] = perm[b];
        perm[b] = temp;
    }
    // Verify
    for (i=1; i<3; i++)
    {
        if (x[perm[i-1]] > x[perm[i]])
        {
#ifdef COMPEARTH_DEBUG_SRC
            fprintf(stderr, "%s: Failed to sort numbers in ascending order\n",
                     __func__);
#endif
            return -1;
        }
    }
    return 0;
}
