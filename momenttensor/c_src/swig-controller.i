%module compearth

%{
#define SWIG_FILE_WITH_INIT
#undef NO_IMPORT_ARRAY
%}

%include "numpy.i"
%include "typemaps.i"
%include "carrays.i"
%init %{
    import_array();
%}

%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int nmt1, int nmt1_2, double * M1),
                                                (int nmt2, int nmt2_2, double * M2),
                                                (int nmt, int nmt_2, double * M)}
/*TODO need to do beter memory management ARGOUTVIEWM should do it but caused segfault */
%apply (int * DIM1, double ** ARGOUTVIEW_ARRAY1) {(int * nomega, double ** omega),
                                                   (int * nM0, double ** M0),
                                                   (int *nMw, double ** Mw),
                                                   (int *nisoPct, double ** isoPct),
                                                   (int *ndevPct, double ** devPct),
                                                   (int *ndcPct, double ** dcPct),
                                                   (int *nclvdPct, double ** clvdPct),
                                                   (int *ngamma, double ** gamma),
                                                   (int *ndelta, double **delta),
                                                   (int *nbeta, double **beta),
                                                   (int *nlambda, double **lambda),
                                                   (int *ntheta, double ** theta)}
%apply (int * DIM1, int * DIM2, double **ARGOUTVIEW_ARRAY2) {(int *nfp1_1, int *nfp1_2, double ** fp1),
                                                              (int *nfp2_1, int *nfp2_2, double ** fp2),
                                                              (int *npAxis_1, int *npAxis_2, double ** pAxis),
                                                              (int *nbAxis_1, int *nbAxis_2, double ** bAxis),
                                                              (int *ntAxis_1, int *ntAxis_2, double ** tAxis),
                                                              (int *nM_1, int *nM_2, double ** M)}
%apply (int * DIM1, int * DIM2, int * DIM3, double **ARGOUTVIEW_ARRAY3) {(int *nU_1, int *nU_2, int *nU_3, double ** U)}
/*XXX %apply (double ARGOUT_ARRAY1[ANY]) {(double M[6]),
                             (double lam[3]),
                             (double U[9])}*/
%apply (int DIM1, double *IN_ARRAY1) {(int nv, double *v),
                                      (int nh, double *h),
                                      (int nw, double *w),
                                      (int nu, double *u),
                                      (int nM0, double * M0),
                                      (int nsigmaIn, double * sigmaIn),
                                      (int ngamma, double * gamma),
                                      (int ndelta, double *delta),
                                      /*TODO (int nbeta, double *beta),*/
                                      (int nkappa, double *kappa),
                                      (int ntheta, double * theta)}

//%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double **omega)}

//XXX%rename (tt2cmt) compearth_tt2cmt;

%inline %{
#include "compearth.h"
#include "compearth_constants.h"

int tt2cmt(int ngamma, double *gamma,
           int ndelta, double *delta,
           int nM0, double *M0,
           int nkappa, double *kappa,
           int ntheta, double *theta,
           int nsigmaIn, double *sigmaIn, 
           int *nM_1, int *nM_2, double ** M,
           int *nlambda,  double ** lambda,
           int *nU_1, int *nU_2, int *nU_3, double ** U)
{
    int ierr=0, i, j ,k;
    double * tmpU=NULL;
    *nM_1 = ngamma;
    *nM_2 = 6;
    *nlambda = ngamma;
    *nU_1 = ngamma;
    *nU_2 = 3;
    *nU_3 = 3;
    *M=(double *)calloc((*nM_1)*(*nM_2), sizeof(double));
    *lambda=(double *)calloc(*nlambda, sizeof(double));
    tmpU=(double *)calloc((*nU_1)*(*nU_2)*(*nU_3), sizeof(double));
    *U=(double *)calloc((*nU_1)*(*nU_2)*(*nU_3), sizeof(double));
    ierr = compearth_TT2CMT(ngamma, gamma, delta, M0, kappa, theta, sigmaIn, *M, *lambda, tmpU);
    /* Transpose matrix */
    for(i = 0; i < *nU_1; i++)
    {
        for(j = 0; j < *nU_2; j++)
        {
            for(k=0; k < j; k++)
            {
                (*U)[i*(*nU_2)*(*nU_3) + j*(*nU_2) + k] = tmpU[i*(*nU_2)*(*nU_3) + k*(*nU_3) + j];
                (*U)[i*(*nU_2)*(*nU_3) + k*(*nU_3) + j] = tmpU[i*(*nU_2)*(*nU_3) + j*(*nU_2) + k];
            }
        }
    }
    free(tmpU);
    return ierr;
}

int rect2lune(int nv, double * v,
              int nw, double * w,
              int *ngamma, double ** gamma,
              int *ndelta, double ** delta)
{
        *ndelta = nw;
        *ngamma = nv;
        *delta=(double *)calloc(*ndelta, sizeof(double));
        *gamma=(double *)calloc(*ngamma, sizeof(double));
        return compearth_rect2lune(nv, v, nw, w, *gamma, *delta);
}
int u2beta(int maxit, int linvType, int nu, double *u, double tol, int *nbeta, double **beta)
{
        *nbeta = nu;
        *beta=(double *)calloc(*nbeta, sizeof(double));
        return compearth_u2beta(nu, maxit,linvType, u, tol, *beta);
}
int v2gamma(int nv, double *v, int *ngamma, double **gamma)
{
   *ngamma = nv;
   *gamma=(double *)calloc(*ngamma, sizeof(double));
   compearth_v2gamma(nv,v, *gamma);
}
int h2theta(int nh, double *h, int *ntheta, double **theta)
{
   *ntheta = nh;
   *theta=(double *)calloc(*ntheta, sizeof(double));
   compearth_h2theta(nh,h, *theta);
}
 int cmt2omega(int nmt1, int nmt1_2, double * M1,
                        int nmt2, int nmt2_2, double * M2,
                        int * nomega, double ** omega )
    {
        if (*omega)
        {
            free(*omega);
        }
        *nomega = nmt2>nmt1?nmt2:nmt1;
        *omega=(double *)calloc(*nomega, sizeof(double));
        /*TODO error if nmt1_2 or nmt2_2 is not 6 */
        return compearth_CMT2omega(nmt1, M1, nmt2, M2, *omega);
    }
  int standardDecomposition( int nmt, int nmt_2, double * M,
                                    enum compearthCoordSystem_enum basis,
                                    int *nM0, double ** M0,
                                    int *nMw, double ** Mw,
                                    int *nfp1_1, int *nfp1_2, double ** fp1,
                                    int *nfp2_1, int *nfp2_2, double ** fp2,
                                    int *npAxis_1, int *npAxis_2, double ** pAxis,
                                    int *nbAxis_1, int *nbAxis_2, double ** bAxis,
                                    int *ntAxis_1, int *ntAxis_2, double ** tAxis,
                                    int *nisoPct, double ** isoPct,
                                    int *ndevPct, double ** devPct,
                                    int *ndcPct, double ** dcPct,
                                    int *nclvdPct, double ** clvdPct)
    {
        /*TODO if nmt_2 != 6 error */
        *nM0 = nmt;
        *M0=(double *)calloc(*nM0, sizeof(double));
        *nMw = nmt;
        *Mw=(double *)calloc(*nMw, sizeof(double));
        *nfp1_1 = nmt;
        *nfp1_2 = 3;
        *fp1=(double *)calloc((*nfp1_1)*(*nfp1_2), sizeof(double));
        *nfp2_1 = nmt;
        *nfp2_2 = 3;
        *fp2=(double *)calloc((*nfp2_1)*(*nfp2_2), sizeof(double));
        *npAxis_1 = nmt;
        *npAxis_2 = 3;
        *pAxis=(double *)calloc((*npAxis_1)*(*npAxis_2), sizeof(double));
        *nbAxis_1 = nmt;
        *nbAxis_2 = 3;
        *bAxis=(double *)calloc((*nbAxis_1)*(*nbAxis_2), sizeof(double));
        *ntAxis_1 = nmt;
        *ntAxis_2 = 3;
        *tAxis=(double *)calloc((*ntAxis_1)*(*ntAxis_2), sizeof(double));
        *nisoPct = nmt;
        *isoPct=(double *)calloc(*nisoPct, sizeof(double));
        *ndevPct = nmt;
        *devPct=(double *)calloc(*ndevPct, sizeof(double));
        *ndcPct = nmt;
        *dcPct=(double *)calloc(*ndcPct, sizeof(double));
        *nclvdPct = nmt;
        *clvdPct=(double *)calloc(*nclvdPct, sizeof(double));
        return  compearth_standardDecomposition(nmt, M, basis, *M0, *Mw, *fp1, *fp2, *pAxis, *bAxis, *tAxis, *isoPct, *devPct, *dcPct, *clvdPct);
    }
%}

//#%ignore (compearth_angleMT)
%include "compearth_constants.h"
