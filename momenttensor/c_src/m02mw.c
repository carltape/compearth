#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Converts from scalar seismic moment to moment magnitude.
 *
 * @param[in] nm      Number of scalar moments.
 * @param[in] imag    KANAMORI_1978 (1) for Kanaomori 1977/1978 formula.
 *                    HARVARD_CMT (2) for GCMT formula.
 * @param[in] M0      Scalar moments in N-m (not dyne-cm).  This is an array
 *                    of dimension [nm].
 *
 * @param[out] Mw     Moment magnitudes.  This is an arrya of dimension [nm].
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape (31-Oct-2007) and translated to C by Ben Baker (2017) 
 *
 * @copyright MIT
 *
 */
int compearth_m02mw(const int nm, const enum magType_enum imag,
                    const double *__restrict__ M0,
                    double *__restrict__ Mw)
{
    const char *fcnm = "compearth_m02mw\0";
    double m0dcm;
    int i;
    const double A = 2.0/(3.0*log(10.0));
    const double K = 0.2*pow(10.0, 16.8);
    const double two_third = 2.0/3.0;
    if (imag == KANAMORI_1978)
    {
        for (i=0; i<nm; i++)
        {
            m0dcm = 1.e7*M0[i];
            Mw[i] = A*log(m0dcm/K);
        }
    }
    else if (imag == HARVARD_CMT)
    {
        for (i=0; i<nm; i++)
        {
            m0dcm = 1.e7*M0[i];
            Mw[i] = two_third*(log10(m0dcm) - 16.1);
        }
    }
    else
    {
        printf("%s: Invalid magnitude type\n", fcnm);
        for (i=0; i<nm; i++){Mw[i] = 0.0;}
    }
    return 0;
} 
