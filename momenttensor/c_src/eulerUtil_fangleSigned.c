#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define COMPEARTH_PRIVATE_DET3X3
#include "compearth.h"

/*!
 * @brief Returns the signed angle (of rotation) between two vectors.
 *
 * @param[in] n       Length of vectors.
 * @param[in] va      Initial vector.  This is an array of dimension [n].
 * @param[in] vb      Rotated vector.  This is an array of dimension [n].
 * @param[in] vnor    "Normal of rotation" vector.  This is an array of
 *                    dimension [n].
 *
 * @param[out] ierr   0 indicates success.
 *
 * @result Signed angle of rotation between va and vb in degrees.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
double compearth_eulerUtil_fangleSigned(const int n,
                                        const double *__restrict__ va,
                                        const double *__restrict__ vb,
                                        const double *__restrict__ vnor,
                                        int *ierr)
{
    double Dmat[9], stheta;
    *ierr = 0;
    stheta = compearth_matlab_fangle(n, va, vb);
    if (fabs(stheta - 180.0) < 0.0)
    {
        stheta = 180.0;
    }
    else
    {
        Dmat[0] = va[0];
        Dmat[1] = va[1];
        Dmat[2] = va[2];
        Dmat[3] = vb[0];
        Dmat[4] = vb[1]; 
        Dmat[5] = vb[2];
        Dmat[6] = vnor[0];
        Dmat[7] = vnor[1];
        Dmat[8] = vnor[2];
        if (det3x3ColumnMajor(Dmat) < 0.0){stheta =-stheta;}
    }
    return stheta;
}
