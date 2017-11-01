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
                                                   (int *ntheta, double ** theta)}
%apply (int * DIM1, int * DIM2, double **ARGOUTVIEW_ARRAY2) {(int *nfp1_1, int *nfp1_2, double ** fp1),
                                                              (int *nfp2_1, int *nfp2_2, double ** fp2),
                                                              (int *npAxis_1, int *npAxis_2, double ** pAxis),
                                                              (int *nbAxis_1, int *nbAxis_2, double ** bAxis),
                                                              (int *ntAxis_1, int *ntAxis_2, double ** tAxis)}
%apply (double ARGOUT_ARRAY1[ANY]) {(double M[6]),
                             (double lam[3]),
                             (double U[9])}
%apply (int DIM1, double *IN_ARRAY1) {(int nv, double *v),
                                      (int nh, double *h),
                                      (int nw, double *w),
                                      (int nu, double *u)}
//%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double **omega)}

%rename (tt2cmt) compearth_tt2cmt;

%inline %{
#include "compearth.h"
#include "compearth_constants.h"

extern int compearth_tt2cmt(const double gamma,
                     const double delta,
                     const double M0,
                     const double kappa,
                     const double theta,
                     const double sigmaIn, 
                     double M[6], double lam[3], double U[9]);

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
extern int compearth_tt2cmt(const double gamma,
                     const double delta,
                     const double M0,
                     const double kappa,
                     const double theta,
                     const double sigmaIn, 
                     double M[6], double lam[3], double U[9]);
%include "compearth_constants.h"
