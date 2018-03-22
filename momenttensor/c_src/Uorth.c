#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define COMPEARTH_PRIVATE_DET3X3 1
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_NORM3 1
#define COMPEARTH_PRIVATE_GEM3 1
#define COMPEARTH_PRIVATE_GEMT3 1
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_lapacke.h>
//#include <mkl_cblas.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <lapacke.h>
//#include <cblas.h>
#endif

#define LWORK 15

/*!
 * @brief Converts a nearly orthogonal matrix into a numerically orthogonal
 *        matrix.
 *
 * @param[in] n        Number of bases.
 * @param[in] itype    Orthogonalization strategy. \n
 *                       CE_ORTH_SVD orthgonalizes with the SVD. \n
 *                       CE_ORTH_TAPE2012 orthogalizes with Tape 2012c 
 *                       Appendix E. \n
 *                       CE_ORTH_QUAT orthogonalizes with quaternions.
 * @param[in] Uin      [3 x 3 x n] set of bases where each [3 x 3] basis is
 *                     in column major order.
 *
 * @param[out] Uout    [3 x 3 x n] set of re-orthgonalized bases where
 *                     each [3 x 3] bais is in column major order.
 * @param[out] dtUin   If NULL then this will be ignored.  \n
 *                     Otherwise, this is an array of dimension [n] holding
 *                     the determinants of the input basis.
 * @param[out] dtUout  If NULL then this will be ignored.  \n
 *                     Otherwise, this is an array of dimension [n] holding
 *                     the determinants of the output basis.
 *                    
 * @result 0 indicates success.
 *
 * @author Carl Tape and converted to C by Ben Baker.
 *
 */
int compearth_Uorth(const int n,
                    const enum ceOrthoType_enum itype,
                    const double *__restrict__ Uin,
                    double *__restrict__ Uout,
                    double *__restrict__ dtUin,
                    double *__restrict__ dtUout)
{
    double Ut[9] __attribute__((aligned(64)));
    double U[9] __attribute__((aligned(64)));
    double Vt[9] __attribute__((aligned(64)));
    double work[LWORK], s[3], p[3], det, normb, normp;
    int i, ierr;
    bool lwantDetIn;
    const double tol = DBL_EPSILON*100.0;
    // Check inputs
    if (n < 1 || Uin == NULL || Uout == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: No matrices\n", __func__);}
        if (Uin == NULL){fprintf(stderr, "%s: Uin is NULL\n", __func__);}
        if (Uout == NULL){fprintf(stderr, "%s: Uout is NULL\n", __func__);}
        return -1;
    }
    lwantDetIn = false;
    if (dtUin != NULL){lwantDetIn = true;}
    // Orthogonalize basis with SVD
    if (itype == CE_ORTH_SVD)
    {
        for (i=0; i<n; i++)
        {
            // This will frequently be called with matrices that may already
            // be sufficiently orthonormal.
            det = det3x3ColumnMajor(&Uin[9*i]);
            if (fabs(det - 1.0) > 1.e-12)
            {
                // Compute SVD
                memcpy(Ut, &Uin[9*i], 9*sizeof(double));
                ierr = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'A', 'A',
                                           3, 3, Ut, 3, s, U, 3,
                                           Vt, 3, work, LWORK);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error computing Ut\n", __func__);
                    return -1;
                }
                // Compute U*V' - note Vt is already computed by SVD
                gemm3_colMajorNoTransNoTrans(U, Vt, &Uout[9*i]);
                // N.B. I could get the determinant here from the singular
                // values but i'll wait
            }
            else
            {
                memcpy(&Uout[9*i], &Uin[9*i], 9*sizeof(double));
            }
            if (lwantDetIn){dtUin[i] = det;}
        }
    }
    // Orthgonalize with suggestion in TapeTape2012c Appendix E
    else if (itype == CE_ORTH_TAPE2012)
    {
        for (i=0; i<n; i++)
        {
            // This will frequently be called with matrices that may already
            // be sufficiently orthonormal.
            det = det3x3ColumnMajor(&Uin[9*i]);
            if (fabs(det - 1.0) > tol)
            {
                cross3(&Uin[9*i], &Uin[9*i+3], p);
                normb = norm3(&Uin[9*i+3]);
                normp = norm3(p);
                Uout[9*i]   = Uin[9*i];
                Uout[9*i+1] = Uin[9*i+1];
                Uout[9*i+2] = Uin[9*i+2];
                Uout[9*i+3] = Uin[9*i+3]/normb;
                Uout[9*i+4] = Uin[9*i+4]/normb;
                Uout[9*i+5] = Uin[9*i+5]/normb;
                Uout[9*i+6] = Uin[9*i+6]/normp;
                Uout[9*i+7] = Uin[9*i+7]/normp;
                Uout[9*i+8] = Uin[9*i+8]/normp;
            }
            else
            {
                memcpy(&Uout[9*i], &Uin[9*i], 9*sizeof(double));
            }
            if (lwantDetIn){dtUin[i] = det;}
        }
    }
    else if (itype == CE_ORTH_QUAT)
    {
        fprintf(stderr, "%s: Error quaternions not yet programmed\n",
                 __func__);
        return -1;
    }
    else if (itype == CE_NO_ORTH)
    {
        memcpy(Uout, Uin, 9*(size_t) n*sizeof(double));
        //cblas_dcopy(9*n, Uin, 1, Uout, 1);
    }
    else
    {
        fprintf(stderr, "%s: Only itype=1 is programmed\n", __func__);
        return -1;
    }
    //if (dtUin != NULL)
    //{
    //    for (i=0; i<n; i++)
    //    {
    //        dtUin[i] = det3x3ColumnMajor(&Uin[9*i]);
    //    }
    //}
    if (dtUout != NULL)
    {
        for (i=0; i<n; i++)
        {
            dtUout[i] = det3x3ColumnMajor(&Uout[9*i]);
        }
    }
    return 0; 
}
