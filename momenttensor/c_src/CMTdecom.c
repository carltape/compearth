#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

static int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ iperm);
static int argsort3(const double *__restrict__ x, int *__restrict__ iperm);

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
 */
int compearth_CMTdecom(const int nmt, const double *__restrict__ M,
                       const int isort,
                       double *__restrict__ lam,
                       double *__restrict__ U)
{
    const char *fcnm = "compearth_CMTdecom\0";
    double Lams[3], LamsAbs[3], Mx[9], Ut[9], work[LWORK], *lamsCopy;
    int perm[3], c, i, ierr, info, r;
    // Error checks
    if (nmt < 1 || M == NULL || lam == NULL || U == NULL)
    {
        if (nmt < 1){printf("%s: Error no mts\n", fcnm);}
        if (M == NULL){printf("%s: Error M is NULL\n", fcnm);}
        if (lam == NULL){printf("%s: Error lam is NULL\n", fcnm);}
        if (U == NULL){printf("%s: Error U is NULL\n", fcnm);}
        return -1;
    }
    if (isort < 1 || isort > 4)
    {
        printf("%s: Error isort=%d is invalid\n", fcnm, isort);
        return -1;
    }
    // Decompose the moment tensors
    for (i=0; i<nmt; i++)
    {
        compearth_Mvec2Mmat(1, M, 1, Mx);
        ierr = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'L', 3, Mx, 3,
                                  Lams, work, LWORK);
        if (ierr != 0)
        {
            printf("%s: Error computing eigenvalues\n", fcnm);
            return -1;
        }
        if (isort == 1)
        {
            argsort3_upDown(Lams, false, perm); 
            lamsCopy = Lams;
        }
        else if (isort == 2)
        {
            argsort3_upDown(Lams, true, perm);
            lamsCopy = Lams; 
        }
        else if (isort == 3 || isort == 4)
        {
            LamsAbs[0] = fabs(Lams[0]);
            LamsAbs[1] = fabs(Lams[1]);
            LamsAbs[2] = fabs(Lams[2]);
            if (isort == 3)
            {
                argsort3_upDown(LamsAbs, false, perm);
            }
            else
            {
                argsort3_upDown(LamsAbs, true, perm);
            }
            lamsCopy = LamsAbs;
        }
        // And permute and save
        for (r=0; r<3; r++)
        {
            lam[3*i+r] = lamsCopy[perm[r]];
        }
        // Permute columns
        for (c=0; c<3; c++)
        {
            for (r=0; r<3; r++)
            {
                Ut[3*c+r] = Mx[3*c+perm[r]];
            }
        }
        // Ensure the determinant is > 0
        ierr = compearth_Udetcheck(1, Ut, &U[9*i]);
        if (ierr != 0)
        {
            printf("%s: Error checking determinan\n", fcnm);
            return -1;
        }
    }
    return 0;
}

static int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ iperm)
{
    int ipermt[3], ierr;
    if (lascend)
    {
        ierr = argsort3(x, ipermt);
        iperm[0] = ipermt[2];
        iperm[1] = ipermt[1];
        iperm[2] = ipermt[0];
    }
    else
    {
        ierr = argsort3(x, iperm); 
    }
    return ierr;
}

static int argsort3(const double *__restrict__ x, int *__restrict__ iperm)
{
    const char *fcnm = "argsort3\0";
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
            printf("%s: Failed to sort numbers in ascending order\n", fcnm);
            return -1;
        }
    }
    return 0;
}
