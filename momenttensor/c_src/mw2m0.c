#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Convert Mw to M0 using published formulas.
 *
 * @param[in] nm      Number of moment magnitudes.
 * @param[in] imag    KANAMORI_1978 (1) for Kanaomori 1977/1978 formula.
 *                    HARVARD_CMT (2) for GCMT formula.
 * @param[in] Mw      Array of moment magnitudes.  This is an array 
 *                    of dimension [nm].
 *
 * @param[out] M0     Corresponding scalar moments in N-m (i.e. not dyne-cm).
 *                    This is an array of dimension [nm].
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker.
 *
 * @copyright MIT
 *
 */
int compearth_mw2m0(const int nm, const enum magType_enum imag,
                    const double *__restrict__ Mw, 
                    double *__restrict__ M0) 
{
    const char *fcnm = "compearth_m02mw\0";
    const double threeOverTwo = 3.0/2.0;
    const double eleven8_m_log10_5em5 =  11.8 - log10(5e-5);
    int i;
    if (imag == KANAMORI_1978)
    {
        for (i=0; i<nm; i++)
        {
            M0[i] = pow(10., threeOverTwo*Mw[i] + eleven8_m_log10_5em5);
        }
    }
    else if (imag == HARVARD_CMT)
    {
        for (i=0; i<nm; i++)
        {
            M0[i] = pow(10.0, threeOverTwo*Mw[i] + 16.1);
        }
    }
    else
    {
        printf("%s: imag must be 1 or 2\n", fcnm);
        return -1;
    }
    // Formulas are designed for units of dyne-cm, not N-m
    for (i=0; i<nm; i++)
    {
        M0[i] = 1.e-7*M0[i];
    }
    return 0;
}
