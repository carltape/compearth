#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/*!
 * @brief Converts x, y, z triples to latitude and longitude.
 *
 * @param[in] n     Number of points.
 * @param[in] x     x locations.  This is an array of dimension [n].
 * @param[in] y     y locations.  This is an array of dimension [n].
 * @param[in] z     z locations.  This is an array of dimension [n].
 * @param[out] lat  Latitudes (degrees).  This is an array of dimension [n].
 * @param[out] lon  Longitudes (degrees).  This is an array of dimesnion [n].
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
void compearth_xyz2latlon(const int n,
                          const double *__restrict__ x,
                          const double *__restrict__ y,
                          const double *__restrict__ z,
                          double *__restrict__ lat,
                          double *__restrict__ lon)
{
    double th[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double ph[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double rho[CE_CHUNKSIZE] __attribute__((aligned(64)));
    const double toDeg = 180.0/M_PI;
    int i, j, nl; 
    for (j=0; j<n; j=j+CE_CHUNKSIZE)
    {
        nl = MIN(CE_CHUNKSIZE, n-j);
        compearth_xyz2tp(nl, &x[j], &y[j], &z[j], th, ph, rho);
        for (i=0; i<nl; i++)
        {
            lat[i] = (M_PI_2 - th[i])*toDeg;
            lon[i] = ph[i]*toDeg;
        }
    }
    return;
}
