%module compearth
%rename (standardDecomposition) compearth_standardDecomposition;

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
%apply (int * DIM1, double ** ARGOUTVIEWM_ARRAY1) {(int * nomega, double ** omega),
                                                   (int * nM0, double ** M0),
                                                   (int *nMw, double ** Mw),
                                                   (int *nisoPct, double ** isoPct),
                                                   (int *ndevPct, double ** devPct),
                                                   (int *ndcPct, double ** dcPct),
                                                   (int *nclvdPct, double ** clvdPct)}
%apply (int * DIM1, int * DIM2, double **ARGOUTVIEWM_ARRAY2) {(int *nfp1_1, int *nfp1_2, double ** fp1),
                                                              (int *nfp2_1, int *nfp2_2, double ** fp2),
                                                              (int *npAxis_1, int *npAxis_2, double ** pAxis),
                                                              (int *nbAxis_1, int *nbAxis_2, double ** bAxis),
                                                              (int *ntAxis_1, int *ntAxis_2, double ** tAxis)}
//%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double **omega)}

%inline %{
#include "compearth.h"
 int cmt2omega(int nmt1, int nmt1_2, double * M1,
                        int nmt2, int nmt2_2, double * M2,
                        int * nomega, double ** omega )
    {
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
