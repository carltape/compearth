#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"

/*!
 * @brief C translation of Carl Tape's utility which converts a moment
 *        tensor, M, from input system defined by i1in to moment tensor,
 *        Mout, to output system defind by i2in
 *
 * @param[in] i1in     coordinate system for M
 *                       = USE (1) -> up, south, east
 *                       = NED (2) -> north, east, down
 *                       = NWU (3) -> north west, up
 *                       = ENU (4) -> east, north, up
 *                       = SEU (5) -> south, east, up
 * @param[in] i2in     coordinate system for Mout
 *                       = USE (1) -> up, south, east
 *                       = NED (2) -> north, east, down
 *                       = NWU (3) -> north west, up
 *                       = ENU (4) -> east, north, up
 *                       = SEU (5) -> south, east, up
 * @param[in] M        input moment tensor in system i1in [6]
 *                     the C indices {0,1,2,3,4,5} correspond to matrix
 *                     indices: {11, 22, 33, 12, 13, 23}
 *
 * @param[out] Mout    corresponding moment tensor now in system i2in [6]
 *                     the C indices {0,1,2,3,4,5} correspond to matrix
 *                     indices: {11, 22, 33, 12, 13, 23}
 * 
 * @result 0 indicates success
 *
 * @date 2016 - Ben Baker converted Carl Tape's convert_MT.m to C
 *
 * @copyright MIT
 *
 */
int compearth_convertMT(const enum compearthCoordSystem_enum i1in,
                        const enum compearthCoordSystem_enum i2in,
                        const double *__restrict__ M,
                        double *__restrict__ Mout)
{
    const char *fcnm = "compearth_convertMT\0";
    int i1, i2;
    Mout[0] = 0.0;
    Mout[1] = 0.0;
    Mout[2] = 0.0;
    Mout[3] = 0.0;
    Mout[4] = 0.0;
    Mout[5] = 0.0;
    // Quick checks
    i1 = (int) i1in;
    i2 = (int) i2in;
    if (i1 < 1 || i1 > 5)
    {
        printf("%s: Error unknown input coordinate system %d\n", fcnm, i1);
        return -1;
    }
    if (i2 < 1 || i2 > 5)
    {
        printf("%s: Error unkonwn output coordinate system %d\n", fcnm, i2);
        return -1;
    }
    // Base case
    if (i1 == i2)
    {
        Mout[0] = M[0];
        Mout[1] = M[1];
        Mout[2] = M[2];
        Mout[3] = M[3];
        Mout[4] = M[4];
        Mout[5] = M[5];
        return 0;
    }
    // Convert
    if (i1 == 1)
    {
        // up-south-east (GCMT) to north-east-down (AkiRichards)
        // (AR, 1980, p. 118)
        if (i2 == 2)
        {
            Mout[0] = M[1]; //  tt -> xx
            Mout[1] = M[2]; //  pp -> yy
            Mout[2] = M[0]; //  rr -> zz
            Mout[3] =-M[5]; // -tp -> xy
            Mout[4] = M[3]; //  rt -> xz
            Mout[5] =-M[4]; // -rp -> yz
       }
       // up-south-east (GCMT) to north-west-up
       else if (i2 == 3)
       {
            Mout[0] = M[1];
            Mout[1] = M[2];
            Mout[2] = M[0];
            Mout[3] = M[5];
            Mout[4] =-M[3];
            Mout[5] =-M[4];
       }
       // up-south-east (GCMT) to east-north-up
       else if (i2 == 4)
       {
            Mout[0] = M[2];
            Mout[1] = M[1];
            Mout[2] = M[0];
            Mout[3] =-M[5];
            Mout[4] = M[4];
            Mout[5] =-M[3];
       }
       // up-south-east (GCMT) to south-east-up
       else if (i2 == 5)
       {
            Mout[0] = M[1];
            Mout[1] = M[2];
            Mout[2] = M[0];
            Mout[3] = M[5];
            Mout[4] = M[3];
            Mout[5] = M[4];
        }
    }
    else if (i1 == 2)
    {
        // north-east-down (AkiRichards) to up-south-east (GCMT) 
        // (AR, 1980, p. 118)
        if (i2 == 1)
        {
            Mout[0] = M[2]; //  zz -> rr
            Mout[1] = M[0]; //  xx -> tt
            Mout[2] = M[1]; //  yy -> pp 
            Mout[3] = M[4]; //  xz -> rt 
            Mout[4] =-M[5]; // -yz -> rp
            Mout[5] =-M[3]; // -xy -> tp
        }
        // north-east-down (AkiRichards) to north-west-up
        else if (i2 == 3)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] =-M[4];
            Mout[5] = M[5];
        }
        // north-east-down (AkiRichards) to east-north-up
        else if (i2 == 4)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] = M[3];
            Mout[4] =-M[5];
            Mout[5] =-M[4];
        }
        // north-east-down (AkiRichards) to south-east-up
        else if (i2 == 5)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] = M[4];
            Mout[5] =-M[5];
        }
    }
    else if (i1 == 3)
    {
        // north-west-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            Mout[0] = M[2];
            Mout[1] = M[0];
            Mout[2] = M[1];
            Mout[3] =-M[4];
            Mout[4] =-M[5];
            Mout[5] = M[3];
        }
        // north-west-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] =-M[4];
            Mout[5] = M[5];
        }
        // north-west-up to east-north-up
        else if (i2 == 4)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] =-M[5];
            Mout[5] = M[4];
        }
        // north-west-up to south-east-up
        else if (i2 == 5)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] = M[3];
            Mout[4] =-M[4];
            Mout[5] =-M[5];
        }
    }
    else if (i1 == 4)
    {
        // east-north-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            Mout[0] = M[2];
            Mout[1] = M[1];
            Mout[2] = M[0];
            Mout[3] =-M[5];
            Mout[4] = M[4];
            Mout[5] =-M[3];
        }
        // east-north-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] = M[3];
            Mout[4] =-M[5];
            Mout[5] =-M[4];
        }
        // east-north-up to north-west-up
        else if (i2 == 3)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] = M[5];
            Mout[5] =-M[4];
        }
        // east-north-up to south-east-up
        else if (i2 == 5)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] =-M[5];
            Mout[5] = M[4];
        }
    }
    else if (i1 == 5)
    {
        // south-east-up to up-south-east (GCMT)
        if (i2 == 1)
        {
            Mout[0] = M[2];
            Mout[1] = M[0];
            Mout[2] = M[1];
            Mout[3] = M[4];
            Mout[4] = M[5];
            Mout[5] = M[3];
        }
        // south-east-up to north-east-down (AkiRichards)
        else if (i2 == 2)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] = M[4];
            Mout[5] =-M[5];
        }
        // south-east-up to north-west-up
        else if (i2 == 3)
        {
            Mout[0] = M[0];
            Mout[1] = M[1];
            Mout[2] = M[2];
            Mout[3] = M[3];
            Mout[4] =-M[4];
            Mout[5] =-M[5];
        }
        // south-east-up to east-north-up
        else if (i2 == 4)
        {
            Mout[0] = M[1];
            Mout[1] = M[0];
            Mout[2] = M[2];
            Mout[3] =-M[3];
            Mout[4] = M[5];
            Mout[5] =-M[4];
        }
    }
    return 0;
}
