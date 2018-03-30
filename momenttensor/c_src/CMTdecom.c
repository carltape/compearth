#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#define COMPEARTH_PRIVATE_UPDOWN_ARGSORT3 1
#define COMPEARTH_PRIVATE_UPDOWN_ABS_ARGSORT3 1
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_lapacke.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <lapacke.h>
#endif

/*
static int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ iperm);
static int argsort3(const double *__restrict__ x, int *__restrict__ iperm);
*/

#define LWORK 102 //(NB+2)*N hence NB=32; must be atleast 3*N-1 which is 8

/*!
 * @brief Decomposes a set of moment tensors into eigenvalues + basis.
 *
 * @param[in] nmt    Number of moment tensors.
 * @param[in] M      6 x n moment tensors with some unspecified basis 
 *                   (e.g., up-south-east) packed
 *                    M = [M11 M22 M33 M12 M13 M23].
 * @param[in] isort  Sorting of eigenvalues. \n
 *                   If 1 then the eigenvalues are sorted highest to lowest. \n
 *                   If 2 then the eigenvalues are sorted lowest to highest. \n
 *                   If 3 then the eigenvalues are sorted by absolute value
 *                   from highest to lowest. \n
 *                   If 4 then the eigenvalues are sorted by absolute value
 *                   from lowest to highest. 
 *
 * @param[out] lam   3 x n set of eigenvalues.
 * @param[out] U     3 x 3 x n set of bases.  The first 3 x 3 matrix is in
 *                   column major order.  This has the same basis as M.
 *  
 * @author Carl Tape.  Converted to C by Ben Baker.
 * 
 * @copyright MIT
 *
 * @bug The
 */
int compearth_CMTdecom(const int nmt, const double *__restrict__ M,
                       const int isort,
                       double *__restrict__ lam,
                       double *__restrict__ U)
{
    double Lams[3], Mx[9], Ut[9], work[LWORK];
    int perm[3], c, i, ierr, info, r;
    // Error checks
    if (nmt < 1 || M == NULL || lam == NULL || U == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: Error no mts\n", __func__);}
        if (M == NULL){fprintf(stderr, "%s: Error M is NULL\n", __func__);}
        if (lam == NULL){fprintf(stderr, "%s: Error lam is NULL\n", __func__);}
        if (U == NULL){fprintf(stderr, "%s: Error U is NULL\n", __func__);}
        return -1;
    }
    if (isort < 1 || isort > 4)
    {
        fprintf(stderr, "%s: Error isort=%d is invalid\n", __func__, isort);
        return -1;
    }
    // Decompose the moment tensors
    for (i=0; i<nmt; i++)
    {
        // Create the symmetric moment tensor matrix from the 6x1 vector
        compearth_Mvec2Mmat(1, M, 1, Mx);
        // Compute eigenvalues in ascending order
        info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'U', 3, Mx, 3,
                                  Lams, work, LWORK);
        if (info != 0)
        {
            fprintf(stderr, "%s: Error computing eigenvalues\n", __func__);
            return -1;
        }
        // Descending
        if (isort == 1)
        {
            argsort3_upDown(Lams, false, perm); 
        }
        else if (isort == 2)
        {
            argsort3_upDown(Lams, true, perm);
        }
        // Descending on absolute value
        else if (isort == 3)
        {
            argsort3_absUpDown(Lams, false, perm);
        }
        else //if (isort == 4)
        {
            argsort3_absUpDown(Lams, true,  perm);
        }
        // And permute and save
        for (r=0; r<3; r++)
        {
            lam[3*i+r] = Lams[perm[r]];
        }
        // Permute columns
        for (c=0; c<3; c++)
        {
            for (r=0; r<3; r++)
            {
                Ut[3*c+r] = Mx[3*perm[c]+r];
            }
        }
        // Ensure the determinant is > 0
        ierr = compearth_Udetcheck(1, Ut, &U[9*i]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error checking determinant\n", __func__);
            return -1;
        }
    }
    return 0;
}

/*
static int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ iperm)
{
    int ipermt[3], ierr;
    if (lascend)
    {
        ierr = argsort3(x, iperm);
    }
    else
    {
        ierr = argsort3(x, ipermt);
        iperm[0] = ipermt[2];
        iperm[1] = ipermt[1];
        iperm[2] = ipermt[0];
    }
    return ierr;
}

static int argsort3(const double *__restrict__ x, int *__restrict__ iperm)
{
    int i, temp;
    const int a = 0; 
    const int b = 1; 
    const int c = 2; 
    // Copy
    iperm[a] = a; 
    iperm[b] = b; 
    iperm[c] = c; 
    if (x[iperm[a]] > x[iperm[c]])
    {
        temp = iperm[c];
        iperm[c] = iperm[a];
        iperm[a] = temp;
    }    
    if (x[iperm[a]] > x[iperm[b]])
    {
        temp = iperm[b];
        iperm[b] = iperm[a];
        iperm[a] = temp;
    }    
    //Now the smallest element is the first one. Just check the 2-nd and 3-rd
    if (x[iperm[b]] > x[iperm[c]])
    {
        temp = iperm[c];
        iperm[c] = iperm[b];
        iperm[b] = temp;
    }    
    // Verify
    for (i=1; i<3; i++) 
    {
        if (x[iperm[i-1]] > x[iperm[i]])
        {
            fprintf(stderr, "%s: Failed to sort numbers in ascending order\n",
                     __func__);
            return -1;
        }
    }
    return 0;
}
*/
