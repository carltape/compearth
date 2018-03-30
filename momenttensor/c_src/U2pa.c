#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define COMPEARTH_PRIVATE_ANTIPODE 1
#define COMPEARTH_PRIVATE_WRAP360 1
#include "compearth.h"

/*!
 * @brief Convert from basis in South, East, Up to plunge/azimuth of
 *        three basis vectors.
 *
 * @param[in] nmt    Number of bases. 
 * @param[in] U      This is a [3 x 3 x nmt] array of bases in SEU coordinates.
 *                   Each [3 x 3] matrix is in column major format.
 * @param[out] pl1   Plunge angles (degrees) corresponding to first
 *                   eigenvector.
 * @param[out] az1   Azimuth (degrees) corresponding to first eigenvector.
 *                   This is an array of dimension [nmt].
 * @param[out] pl2   Plunge (degrees) corresponding to second eigenvector.
 *                   This is an array of dimension [nmt]
 * @param[out] az2   Azimuth (degrees) corresponding to second eigenvector.
 *                   This is an array of dimension [nmt].
 * @param[out] pl3   Plunge (degrees) corresponding to third eigenvector.
 *                   This is an array of dimension [nmt].
 * @param[out] az3   Azimuth (degrees) corresponding to third  eigenvector.
 *                   This is an array of dimension [nmt].
 *
 * @result 0 indicates success.
 * 
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_U2pa(const int nmt, const double *__restrict__ U,
                   double *__restrict__ pl1,
                   double *__restrict__ az1,
                   double *__restrict__ pl2,
                   double *__restrict__ az2,
                   double *__restrict__ pl3,
                   double *__restrict__ az3)
{
    double Uin[9], Ut[9],
           lat, lon, lat1, lon1, lat2, lon2, lat3, lon3;
    int i, ierr;
    ierr = 0; 
    for (i=0; i<nmt; i++)
    {
        //memcpy(Ut, &U[9*i], 9*sizeof(double)); 
        ierr = compearth_Uorth(1, CE_ORTH_SVD, &U[9*i], Ut, NULL, NULL);
        ierr = compearth_Udetcheck(1, Ut, Uin);
        // Extract column vectors and convert to lat/lon
        //printf("p1: %f %f %f\n", Uin[0], Uin[1], Uin[2]);
        //printf("p2: %f %f %f\n", Uin[3], Uin[4], Uin[5]);
        //printf("p3: %f %f %f\n", Uin[6], Uin[7], Uin[8]);
        compearth_xyz2latlon(1, &Uin[0], &Uin[1], &Uin[2], &lat1, &lon1);
        compearth_xyz2latlon(1, &Uin[3], &Uin[4], &Uin[5], &lat2, &lon2);
        compearth_xyz2latlon(1, &Uin[6], &Uin[7], &Uin[8], &lat3, &lon3);
        //printf("lat1/lon1 = %f %f\n", lat1,lon1);
        //printf("lat2/lon2 = %f %f\n", lat2,lon2);
        //printf("lat3/lon3 = %f %f\n", lat3,lon3);
        // plunge must be 0 to 90
        if (Uin[2] < 0.0) // p1(3)
        {
            lat = lat1;
            lon = lon1; 
            antipode(lat, lon, &lat1, &lon1, true);
        }
        if (Uin[5] < 0.0) // p2(3)
        {
            lat = lat2;
            lon = lon2; 
            antipode(lat, lon, &lat2, &lon2, true);
        }
        if (Uin[8] < 0.0) // p1(3)
        {
            lat = lat3;
            lon = lon3;
            antipode(lat, lon, &lat3, &lon3, true);
        }
        //printf("%f %f\n", lat1, lon1);
        pl1[i] = lat1;
        pl2[i] = lat2;
        pl3[i] = lat3;
        az1[i] = wrap360(-lon1);
        az2[i] = wrap360(-lon2);
        az3[i] = wrap360(-lon3);
    }
    return ierr;
}
