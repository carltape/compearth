#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compearth.h"
/*
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
*/

/*!
 * @brief Converts a moment tensor, M, from input system defined by i1in
 *        to moment tensor, Mout, to output system defind by i2in.
 *
 * @param[in] nmt      Number of moment tensors.
 * @param[in] i1in     Coordinate system for M: \n
 *                       = CE_USE (1) -> up, south, east \n
 *                       = CE_NED (2) -> north, east, down \n
 *                       = CE_NWU (3) -> north west, up \n
 *                       = CE_ENU (4) -> east, north, up \n
 *                       = CE_SEU (5) -> south, east, up \n
 * @param[in] i2in     Coordinate system for Mout: \n
 *                       = CE_USE (1) -> up, south, east \n
 *                       = CE_NED (2) -> north, east, down \n
 *                       = CE_NWU (3) -> north west, up \n
 *                       = CE_ENU (4) -> east, north, up \n
 *                       = CE_SEU (5) -> south, east, up \n
 * @param[in] M        Input moment tensor in system i1in.  This is an
 *                     an array of dimension [6 x nmt].
 *                     The C indices {0,1,2,3,4,5} correspond to matrix
 *                     indices: {11, 22, 33, 12, 13, 23}.
 *
 * @param[out] Mout    Corresponding moment tensor now in system i2in.
 *                     This is an array of dimension [6 x nmt].
 *                     the C indices {0,1,2,3,4,5} correspond to matrix
 *                     indices: {11, 22, 33, 12, 13, 23}
 * 
 * @result 0 indicates success.
 *
 * @date 2016 - Ben Baker converted Carl Tape's convert_MT.m to C
 *
 * @copyright MIT
 *
 */
int compearth_convertMT(const int nmt,
                        const enum compearthCoordSystem_enum i1in,
                        const enum compearthCoordSystem_enum i2in,
                        const double *__restrict__ M,
                        double *__restrict__ Mout)
{
    int i, i1, i2;
    // Check the inputs to avoid seg faults
    if (nmt < 1 || M == NULL || Mout == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (M == NULL){fprintf(stderr, "%s: Error M is NULL\n", __func__);}
        if (Mout == NULL)
        {
            fprintf(stderr, "%s: Error Mout is NULL\n", __func__);
        }
        return -1;
    }
    for (i=0; i<6*nmt; i++){Mout[i] = 0.0;}
    // Quick checks
    i1 = (int) i1in;
    i2 = (int) i2in;
    if (i1 < 1 || i1 > 5)
    {
        fprintf(stderr, "%s: Error unknown input coordinate system %d\n",
                __func__, i1);
        return -1;
    }
    if (i2 < 1 || i2 > 5)
    {
        fprintf(stderr, "%s: Error unkonwn output coordinate system %d\n",
                __func__, i2);
        return -1;
    }
    // Base case
    if (i1 == i2)
    {
        memcpy(Mout, M, (size_t) nmt*6*sizeof(double));
        //cblas_dcopy(6*nmt, M, 1, Mout, 1);
        return 0;
    }
    // Convert
    if (i1 == 1)
    {
        // up-south-east (GCMT) to north-east-down (AkiRichards)
        // (AR, 1980, p. 118)
        if (i2 == 2)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1]; //  tt -> xx
                Mout[6*i+1] = M[6*i+2]; //  pp -> yy
                Mout[6*i+2] = M[6*i+0]; //  rr -> zz
                Mout[6*i+3] =-M[6*i+5]; // -tp -> xy
                Mout[6*i+4] = M[6*i+3]; //  rt -> xz
                Mout[6*i+5] =-M[6*i+4]; // -rp -> yz
            }
       }
       // up-south-east (GCMT) to north-west-up
       else if (i2 == 3)
       {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+2];
                Mout[6*i+2] = M[6*i+0];
                Mout[6*i+3] = M[6*i+5];
                Mout[6*i+4] =-M[6*i+3];
                Mout[6*i+5] =-M[6*i+4];
            }
       }
       // up-south-east (GCMT) to east-north-up
       else if (i2 == 4)
       {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+2];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+0];
                Mout[6*i+3] =-M[6*i+5];
                Mout[6*i+4] = M[6*i+4];
                Mout[6*i+5] =-M[6*i+3];
            }
       }
       // up-south-east (GCMT) to south-east-up
       else if (i2 == 5)
       {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+2];
                Mout[6*i+2] = M[6*i+0];
                Mout[6*i+3] = M[6*i+5];
                Mout[6*i+4] = M[6*i+3];
                Mout[6*i+5] = M[6*i+4];
            }
        }
    }
    else if (i1 == 2)
    {
        // north-east-down (AkiRichards) to up-south-east (GCMT) 
        // (AR, 1980, p. 118)
        if (i2 == 1)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+2]; //  zz -> rr
                Mout[6*i+1] = M[6*i+0]; //  xx -> tt
                Mout[6*i+2] = M[6*i+1]; //  yy -> pp 
                Mout[6*i+3] = M[6*i+4]; //  xz -> rt 
                Mout[6*i+4] =-M[6*i+5]; // -yz -> rp
                Mout[6*i+5] =-M[6*i+3]; // -xy -> tp
            }
        }
        // north-east-down (AkiRichards) to north-west-up
        else if (i2 == 3)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] =-M[6*i+4];
                Mout[6*i+5] = M[6*i+5];
            }
        }
        // north-east-down (AkiRichards) to east-north-up
        else if (i2 == 4)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] = M[6*i+3];
                Mout[6*i+4] =-M[6*i+5];
                Mout[6*i+5] =-M[6*i+4];
            }
        }
        // north-east-down (AkiRichards) to south-east-up
        else if (i2 == 5)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] = M[6*i+4];
                Mout[6*i+5] =-M[6*i+5];
            }
        }
    }
    else if (i1 == 3)
    {
        // north-west-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+2];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+1];
                Mout[6*i+3] =-M[6*i+4];
                Mout[6*i+4] =-M[6*i+5];
                Mout[6*i+5] = M[6*i+3];
            }
        }
        // north-west-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] =-M[6*i+4];
                Mout[6*i+5] = M[6*i+5];
            }
        }
        // north-west-up to east-north-up
        else if (i2 == 4)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] =-M[6*i+5];
                Mout[6*i+5] = M[6*i+4];
            }
        }
        // north-west-up to south-east-up
        else if (i2 == 5)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] = M[6*i+3];
                Mout[6*i+4] =-M[6*i+4];
                Mout[6*i+5] =-M[6*i+5];
            }
        }
    }
    else if (i1 == 4)
    {
        // east-north-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+2];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+0];
                Mout[6*i+3] =-M[6*i+5];
                Mout[6*i+4] = M[6*i+4];
                Mout[6*i+5] =-M[6*i+3];
            }
        }
        // east-north-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] = M[6*i+3];
                Mout[6*i+4] =-M[6*i+5];
                Mout[6*i+5] =-M[6*i+4];
            }
        }
        // east-north-up to north-west-up
        else if (i2 == 3)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] = M[6*i+5];
                Mout[6*i+5] =-M[6*i+4];
            }
        }
        // east-north-up to south-east-up
        else if (i2 == 5)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] =-M[6*i+5];
                Mout[6*i+5] = M[6*i+4];
            }
        }
    }
    else if (i1 == 5)
    {
        // south-east-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+2];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+1];
                Mout[6*i+3] = M[6*i+4];
                Mout[6*i+4] = M[6*i+5];
                Mout[6*i+5] = M[6*i+3];
            }
        }
        // south-east-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] = M[6*i+4];
                Mout[6*i+5] =-M[6*i+5];
            }
        }
        // south-east-up to north-west-up
        else if (i2 == 3)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+0];
                Mout[6*i+1] = M[6*i+1];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] = M[6*i+3];
                Mout[6*i+4] =-M[6*i+4];
                Mout[6*i+5] =-M[6*i+5];
            }
        }
        // south-east-up to east-north-up
        else if (i2 == 4)
        {
            for (i=0; i<nmt; i++)
            {
                Mout[6*i+0] = M[6*i+1];
                Mout[6*i+1] = M[6*i+0];
                Mout[6*i+2] = M[6*i+2];
                Mout[6*i+3] =-M[6*i+3];
                Mout[6*i+4] = M[6*i+5];
                Mout[6*i+5] =-M[6*i+4];
            }
        }
    }
    return 0;
}
