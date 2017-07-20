#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes rotation matrix for given axis and angle
 *
 * @param[in] n        number of matrices
 * @param[in] xdeg     rotation angles (degrees) [n]
 * @param[in] ixyz     if 1 then rotate about x.
 *                     if 2 then rotate about y.
 *                     if 3 then rotate about z.
 *
 * @param[out] R       rotation matrix in column major format [3 x 3 x n]
 *
 * @result 0 indicates success
 *
 * @date 2016 - Ben Baker converted Carl Tape's rotmat.m to C
 *
 * @copyright MIT
 *
 */
int compearth_eulerUtil_rotmat(const int n, const double *__restrict__ xdeg,
                               const int ixyz, double *__restrict__ R)
{
    const char *fcnm = "compearth_eulerUtil_rotmat\0";
    double cosx, sinx;
    const double pi180 = M_PI/180.0;
    int i;
    if (ixyz == 1)
    {
        for (i=0; i<n; i++)
        {
            cosx = cos(xdeg[i]*pi180);
            sinx = sin(xdeg[i]*pi180);
            R[9*i+0] = 1.0;
            R[9*i+1] = 0.0;
            R[9*i+2] = 0.0;
            R[9*i+3] = 0.0;
            R[9*i+4] = cosx;
            R[9*i+5] = sinx;
            R[9*i+6] = 0.0; 
            R[9*i+7] =-sinx;
            R[9*i+8] = cosx;
        }
    }
    else if (ixyz == 2)
    {
        for (i=0; i<n; i++)
        {
            cosx = cos(xdeg[i]*pi180);
            sinx = sin(xdeg[i]*pi180);
            R[9*i+0] = cosx;
            R[9*i+1] = 0.0;
            R[9*i+2] =-sinx;
            R[9*i+3] = 0.0;
            R[9*i+4] = 1.0;
            R[9*i+5] = 0.0;
            R[9*i+6] = sinx;
            R[9*i+7] = 0.0;
            R[9*i+8] = cosx;
        }
    }
    else if (ixyz == 3)
    {
        for (i=0; i<n; i++)
        {
            cosx = cos(xdeg[i]*pi180);
            sinx = sin(xdeg[i]*pi180);
            R[9*i+0] = cosx;
            R[9*i+1] = sinx;
            R[9*i+2] = 0.0;
            R[9*i+3] =-sinx;
            R[9*i+4] = cosx;
            R[9*i+5] = 0.0;
            R[9*i+6] = 0.0;
            R[9*i+7] = 0.0;
            R[9*i+8] = 1.0;
        }
    }
    else
    {
        printf("%s: Invalid rotation %d; only 1, 2, or 3\n", fcnm, ixyz);
        return -1;
    }    
    return 0;
}
