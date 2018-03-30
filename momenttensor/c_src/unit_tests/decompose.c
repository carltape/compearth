#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"
#include "cmopad.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_cblas.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <cblas.h>
#endif

struct mt_struct
{
    double mt[6];
    double tax[3];
    double bax[3];
    double pax[3];
    double sdr1[3]; //strike dip rake; fault plane 2
    double sdr2[3]; //strike dip rake; fault plane 1
    double mag;
    double exp;
};

/*
int compearth_standardDecomposition(const int nmt,
                                    const double *__restrict__ M,
                                    enum compearthCoordSystem_enum basis,
                                    double *__restrict__ M0,
                                    double *__restrict__ Mw,
                                    double *__restrict__ fp1,
                                    double *__restrict__ fp2,
                                    double *__restrict__ pAxis,
                                    double *__restrict__ bAxis,
                                    double *__restrict__ tAxis,
                                    double *__restrict__ isoPct,
                                    double *__restrict__ devPct,
                                    double *__restrict__ dcPct,
                                    double *__restrict__ clvdPct);
*/
int decompose(const int nmt, const double mtUSE[6]);

int decompose(const int nmt, const double mtUSE[6])
{
    //double K[3], N[3], S[3], lam[3], U[9], Uuse[9], b[3], p[3], t[3], 
    //double gamma, delta, M0, kappa, theta, sigma, thetadc,
    //       k2, d2, s2;
    struct cmopad_struct src;
    const char *fcnm = "decompose\0";
    enum cmopad_basis_enum cin, cloc;
    double M[3][3], pax[3], bax[3], tax[3];
    int i, ierr;
    const int iverb = 0;
/*
    //------------------------------------------------------------------------//
    //                       Do the cMoPaD decomposition                      //
    //------------------------------------------------------------------------//
    // Copy moment tensor
    M[0][0] = mtUSE[0];           //Mrr
    M[1][1] = mtUSE[1];           //Mtt
    M[2][2] = mtUSE[2];           //Mpp
    M[0][1] = M[1][0] = mtUSE[3]; //Mrt
    M[0][2] = M[2][0] = mtUSE[4]; //Mrp
    M[1][2] = M[2][1] = mtUSE[5]; //Mtp
    // Mopad works in North, East, Down frame but inversion is Up, South, East
    cin = USE;
    cloc = NED;
    cmopad_basis_transformMatrixM33(M, cin, cloc); //USE -> NED
    // Compute the isotropic, CLVD, DC decomposition 
    ierr = cmopad_standardDecomposition(M, &src);
    if (ierr != 0)
    {
        printf("%s: Error in decomposition!\n", fcnm);
        return EXIT_FAILURE;
    }
    // Compute the princple axes with corresponding strikes, dips, and rakes 
    ierr = cmopad_MT2PrincipalAxisSystem(iverb, &src);
    if (ierr != 0)
    {
        printf("%s: Error computing principal axis\n",fcnm);
        return EXIT_FAILURE;
    }
    // Compute the pressure, null, and, tension principal axes as len,az,plunge
    cin = NED; //MoPaD is in north, east, down 
    ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[0],
                                            src.p_axis,    pax);
    if (ierr != 0)
    {
        printf("%s: Error converting pax\n",fcnm);
        return EXIT_FAILURE;
    }
    ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[1],
                                            src.null_axis, bax);
    if (ierr != 0){
        printf("%s: Error converting bax\n",fcnm);
        return EXIT_FAILURE;
    }
    ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[2],
                                            src.t_axis,    tax);
    if (ierr != 0)
    {
        printf("%s: Error converting tax\n",fcnm);
        return EXIT_FAILURE;
    }
*/
    //printf("%f %f %f %f %f %f\n",
    //       mtUSE[0], mtUSE[1], mtUSE[2], mtUSE[3], mtUSE[4], mtUSE[5]);
    //ierr = compearth_CMT2TT(1, mtUSE, false,
    //                        &gamma, &delta, &M0,
    //                        &kappa, &theta, &sigma,
    //                        K, N, S, &thetadc,  
    //                        lam, U);
    //if (ierr != 0)
    //{
    //    fprintf(stderr, "%s: CMT2TT failed\n", __func__);
    //    return EXIT_FAILURE;
    //}
    /*
    // Need to switch from SEU to USE.  This requires computing
    // R U inv(R) where:
    //  R = [0 0 1
    //       1 0 0
    //       0 1 0]
    // Carrying through the multiplication yields the permutation
    //  U_{use} = [U_{33} U_{31} U_{32}
    //             U_{13} U_{11} U_{12}
    //             U_{23} U_{21} U_{22}]
    //
    // Permute the column major matrix
    Uuse[0] = U[8];
    Uuse[1] = U[6];
    Uuse[2] = U[7];

    Uuse[3] = U[2];
    Uuse[4] = U[0];
    Uuse[5] = U[1];

    Uuse[6] = U[5];
    Uuse[7] = U[3];
    Uuse[8] = U[4];
printf("%f %f %f\n%f %f %f\n%f %f %f\n",
       U[0], U[3], U[6],
       U[1], U[4], U[7],
       U[2], U[5], U[8]);
*/
    //ierr = compearth_U2pa(1, U,
    //                      &p[2], &p[1], &b[2], &b[1], &t[2], &t[1]);
    //if (ierr != 0)
    //{
    //    fprintf(stderr, "%s: Error converting u to plunge/azimuth\n", 
    //            __func__);
    //    return EXIT_FAILURE;
    //}
    //// Fill in the eigenvalues
    //p[0] = lam[0];
    //b[0] = lam[1];
    //t[0] = lam[2];
    //printf("p: %e %f %f\n", p[0], p[1], p[2]);
    //printf("b: %e %f %f\n", b[0], b[1], b[2]);
    //printf("t: %e %f %f\n", t[0] ,t[1], t[2]); 
    // Compute the auxiliary fault plane
    //ierr = compearth_auxiliaryPlane(1, &kappa, &theta, &sigma,
    //                                &k2, &d2, &s2);
    // Just define some convention where smaller strike goes first
    //printf("%e %f %f %f\n", M0, kappa, theta, sigma);
    //printf("%f %f %f\n", k2, d2, s2);
    //------------------------------------------------------------------------//
    //                         do the standard decomposition                  //
    //------------------------------------------------------------------------//
    //double pAxis[3], bAxis[3], tAxis[3], fp1[3], fp2[3], isoPct, dcPct, devPct, clvdPct, Mw;
    double *pAxis, *bAxis, *tAxis, *fp1, *fp2, *isoPct, *dcPct, *devPct, *clvdPct, *Mw, *M0;
    M0 = (double *) calloc((size_t) nmt, sizeof(double));
    Mw = (double *) calloc((size_t) nmt, sizeof(double));
    fp1 = (double *) calloc((size_t) (3*nmt), sizeof(double));
    fp2 = (double *) calloc((size_t) (3*nmt), sizeof(double));
    pAxis = (double *) calloc((size_t) (3*nmt), sizeof(double));
    bAxis = (double *) calloc((size_t) (3*nmt), sizeof(double));
    tAxis = (double *) calloc((size_t) (3*nmt), sizeof(double));
    isoPct = (double *) calloc((size_t) nmt, sizeof(double));
    devPct = (double *) calloc((size_t) nmt, sizeof(double));
    dcPct = (double *) calloc((size_t) nmt, sizeof(double));
    clvdPct = (double *) calloc((size_t) nmt, sizeof(double));
    ierr = compearth_standardDecomposition(nmt, mtUSE, CE_USE, 
                                           M0, Mw, fp1, fp2,
                                           pAxis, bAxis, tAxis,
                                           isoPct, devPct, dcPct, clvdPct);
    if (ierr != 0)
    {
        fprintf(stderr, "Failed to compute standard decomposition\n");
        return EXIT_FAILURE;
    }
    for (i=0; i<nmt; i++)
    {
        //--------------------------------------------------------------------//
        //                       Do the cMoPaD decomposition                  //
        //--------------------------------------------------------------------//
        // Copy moment tensor
        M[0][0] = mtUSE[6*i+0];           //Mrr
        M[1][1] = mtUSE[6*i+1];           //Mtt
        M[2][2] = mtUSE[6*i+2];           //Mpp
        M[0][1] = M[1][0] = mtUSE[6*i+3]; //Mrt
        M[0][2] = M[2][0] = mtUSE[6*i+4]; //Mrp
        M[1][2] = M[2][1] = mtUSE[6*i+5]; //Mtp
        // Mopad works in North, East, Down
        cin = USE;
        cloc = NED;
        cmopad_basis_transformMatrixM33(M, cin, cloc); //USE -> NED
        // Compute the isotropic, CLVD, DC decomposition 
        ierr = cmopad_standardDecomposition(M, &src);
        if (ierr != 0)
        {   
            printf("%s: Error in decomposition!\n", fcnm);
            return EXIT_FAILURE; 
        }   
        // Compute the princple axes and strikes, dips, and rakes 
        ierr = cmopad_MT2PrincipalAxisSystem(iverb, &src);
        if (ierr != 0)
        {   
            printf("%s: Error computing principal axis\n",fcnm);
            return EXIT_FAILURE; 
        }   
        // Compute the pressure, null, and, tension principal axes as
        // len,az,plunge
        cin = NED; //MoPaD is in north, east, down 
        ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[0],
                                                src.p_axis,    pax);
        if (ierr != 0)
        {   
            printf("%s: Error converting pax\n",fcnm);
            return EXIT_FAILURE; 
        }   
        ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[1],
                                                src.null_axis, bax);
        if (ierr != 0)
        { 
            printf("%s: Error converting bax\n",fcnm);
            return EXIT_FAILURE; 
        }   
        ierr = cmopad_Eigenvector2PrincipalAxis(cin, src.eig_pnt[2],
                                                src.t_axis,    tax);
        if (ierr != 0)
        {   
            printf("%s: Error converting tax\n",fcnm);
            return EXIT_FAILURE; 
        }
        // Check it
        //printf("%d %d %e\n", i, i%9, Mw[i]);
        if (fabs(M0[i] - src.seismic_moment)/src.seismic_moment > 1.e-10)
        {
            fprintf(stderr, "failed M0 %f %f %e %d\n", M0[i], src.seismic_moment,
                    fabs(M0[i] - src.seismic_moment)/src.seismic_moment, i);
            return EXIT_FAILURE;
        }
        if (fabs(Mw[i] - src.moment_magnitude) > 1.e-10)
        {
            fprintf(stderr, "failed Mw %f %f %d\n", Mw[i], src.moment_magnitude, i);
            return EXIT_FAILURE;
        }
        if (fabs(isoPct[i] - src.ISO_percentage) > 1.e-10)
        {   
            fprintf(stderr, "Failed iso pct %f %f\n",
                    isoPct[i], src.ISO_percentage);
            return EXIT_FAILURE;
        }
        if (fabs(dcPct[i] - src.DC_percentage) > 1.e-10)
        {
            fprintf(stderr, "Failed DC pct %f %f %d\n",
                    dcPct[i], src.DC_percentage, i);
            return EXIT_FAILURE;
        }
        if (fabs(devPct[i] - src.DEV_percentage) > 1.e-10)
        {
            fprintf(stderr, "Failed dev pct %f %f\n",
                    devPct[i], src.DEV_percentage);
            return EXIT_FAILURE;
        }
        if (fabs(clvdPct[i] - src.CLVD_percentage) > 1.e-10)
        {
            fprintf(stderr, "Failed clvd pct %f %f\n",
                    devPct[i], src.CLVD_percentage);
            return EXIT_FAILURE;
        }
        // Swap fault planes
        if (fabs(fp1[3*i+0] - src.fp1[0]) > 1.e-10)
        {
            if (fabs(fp1[3*i+0] - src.fp2[0]) > 1.e-10)
            {
                fprintf(stderr, "Failed strike 1 %f %f\n", fp1[3*i+0], src.fp2[0]);
                return EXIT_FAILURE;
            }
            if (fabs(fp1[3*i+1] - src.fp2[1]) > 1.e-10)
            {
                fprintf(stderr, "Failed dip 1 %f %f\n", fp1[3*i+1], src.fp2[1]);
                return EXIT_FAILURE;
            }
            if (fabs(fp1[3*i+2] - src.fp2[2]) > 1.e-10)
            {
               fprintf(stderr, "Failed rake 1 %f %f\n", fp1[3*i+2], src.fp2[2]);
               return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+0] - src.fp1[0]) > 1.e-10)
            {
                fprintf(stderr, "Failed strike 2 %f %f\n", fp2[3*i+0], src.fp1[0]);
                return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+1] - src.fp1[1]) > 1.e-10)
            {
                fprintf(stderr, "Failed dip 2 %f %f\n", fp2[3*i+1], src.fp1[1]);
                return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+2] - src.fp1[2]) > 1.e-10)
            {
                fprintf(stderr, "Failed rake 2 %f %f\n", fp2[3*i+2], src.fp1[2]);
                return EXIT_FAILURE;
            }
        }
        else
        {
            if (fabs(fp1[3*i+0] - src.fp1[0]) > 1.e-10)
            {
                fprintf(stderr, "Failed strike 1 %f %f\n",
                        fp1[3*i+0], src.fp1[0]);
                return EXIT_FAILURE;
            }
            if (fabs(fp1[3*i+1] - src.fp1[1]) > 1.e-10)
            {
                fprintf(stderr, "Failed dip 1 %f %f\n",
                        fp1[3*i+1], src.fp1[1]);
                return EXIT_FAILURE;
            }
            if (fabs(fp1[3*i+2] - src.fp1[2]) > 1.e-10)
            {
                fprintf(stderr, "Failed rake 1 %f %f\n",
                        fp1[3*i+2], src.fp1[2]);
                return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+0] - src.fp2[0]) > 1.e-10)
            {
                fprintf(stderr, "Failed strike 2 %f %f\n",
                       fp2[3*i+0], src.fp2[0]);
                return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+1] - src.fp2[1]) > 1.e-10)
            {
                fprintf(stderr, "Failed dip 2 %f %f\n",
                        fp2[3*i+1], src.fp2[1]);
                return EXIT_FAILURE;
            }
            if (fabs(fp2[3*i+2] - src.fp2[2]) > 1.e-10)
            {
                fprintf(stderr, "Failed rake 2 %f %f\n",
                        fp2[3*i+2], src.fp2[2]);
                return EXIT_FAILURE;
            }
        }
        if (fabs(pAxis[3*i+0] - pax[0]) > 1.e-10 ||
            fabs(pAxis[3*i+1] - pax[1]) > 1.e-10 ||
            fabs((pAxis[3*i+2] - pax[2])/pax[2]) > 1.e-10)
        {
            fprintf(stderr, "failed pressure axix: %f %f %f %f %f %f\n",
                    pAxis[0], pAxis[1], pAxis[2], pax[0], pax[1], pax[2]);
            return EXIT_FAILURE;
        }
        if (fabs(bAxis[3*i+0] - bax[0]) > 1.e-10 ||
            fabs(bAxis[3*i+1] - bax[1]) > 1.e-10 ||
            fabs((bAxis[3*i+2] - bax[2])/bax[2]) > 1.e-10)
        {
            fprintf(stderr, "Failed null axis: %f %f %f %f %f %f\n",
                    bAxis[0], bAxis[1], bAxis[2], bax[0], bax[1], bax[2]);
            return EXIT_FAILURE;
        }
        if (fabs(tAxis[3*i+0] - tax[0]) > 1.e-10 ||
            fabs(tAxis[3*i+1] - tax[1]) > 1.e-10 ||
            fabs((tAxis[3*i+2] - tax[2])/tax[2]) > 1.e-10)
        {
            fprintf(stderr, "Failed tension axis %f %f %f %f %f %f\n",
                    tAxis[0], tAxis[1], tAxis[2], tax[0], tax[1], tax[2]);
            return EXIT_FAILURE;
        }
        if (nmt == 1)
        {
            printf("Summary:\n");
            printf("Mw: %f\n", Mw[i]);
            printf("(strike,dip,rake): (%f,%f,%f)\n",
                   fp1[3*i+0], fp1[3*i+1], fp1[3*i+2]);
            printf("(strike,dip,rake): (%f,%f,%f)\n",
                   fp2[3*i+0], fp2[3*i+1], fp2[3*i+2]);
            printf("iso pct: %f\n", isoPct[i]);
            printf("dc pct: %f\n", dcPct[i]);
            printf("dev pct: %f\n", devPct[i]); 
            printf("clvd pct: %f\n", clvdPct[i]);
            printf("\n");
        }
    }
    free(M0);
    free(Mw);
    free(fp1);
    free(fp2);
    free(pAxis);
    free(bAxis);
    free(tAxis);
    free(isoPct);
    free(devPct);
    free(dcPct);
    free(clvdPct);
    return EXIT_SUCCESS;
} 

int main()
{
    struct mt_struct mt;
    const int nmt = 21*CE_CHUNKSIZE - 1;
    double *mtTest = (double *) calloc((size_t) (6*nmt), sizeof(double));
             //mtTest[2*6*CE_CHUNKSIZE];
    double xscal;
    int ierr, jmt;
    const int n6 = 6;
    const int incx = 1;
    //-----------------------------------1------------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 5.1;
    mt.mt[0] =-4.360;
    mt.mt[1] = 3.200;
    mt.mt[2] = 1.170;
    mt.mt[3] =-0.794;
    mt.mt[4] = 1.650;
    mt.mt[5] = 2.060;
    mt.sdr1[0] = 77.0;
    mt.sdr1[1] = 47.0;
    mt.sdr1[2] =-62.0;
    mt.sdr2[0] = 219.0;
    mt.sdr2[1] = 50.0;
    mt.sdr2[2] =-117.0;
    mt.tax[0] = 4.49;
    mt.tax[1] = 1.0;
    mt.tax[2] = 328.0;
    mt.bax[0] = 0.56;
    mt.bax[1] = 20.0;
    mt.bax[2] = 237.0;
    mt.pax[0] =-5.04;
    mt.pax[1] = 70.0;
    mt.pax[2] = 61.0;
    xscal = pow(10.0, mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[0], 1);
    // Perform decomposition
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposiing 1\n");
        return EXIT_FAILURE;
    }
    //--------------------------------------2----------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 5.0;
    mt.mt[0] =-4.480;
    mt.mt[1] = 1.440; 
    mt.mt[2] = 3.040;
    mt.mt[3] = 0.603; 
    mt.mt[4] = 0.722;
    mt.mt[5] =-1.870;
    mt.sdr1[0] =316.0; 
    mt.sdr1[1] = 44.0;
    mt.sdr1[2] =-105.0; 
    mt.sdr2[0] = 157.0;
    mt.sdr2[1] = 48.0;
    mt.sdr2[2] =-76.0;
    mt.tax[0] = 4.28;
    mt.tax[1] = 2.0;
    mt.tax[2] = 237.0;
    mt.bax[0] = 0.37;
    mt.bax[1] = 10.0;
    mt.bax[2] = 327.0;
    mt.pax[0] =-4.66;
    mt.pax[1] = 79.0;
    mt.pax[2] =137.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[6], 1);
    // Perform decomposition
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {   
        fprintf(stderr, "Error decomposiing 1\n");
        return EXIT_FAILURE;
    }
    //------------------------------------3-----------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 4.9;
    mt.mt[0] =-2.460;
    mt.mt[1] = 0.207; 
    mt.mt[2] = 2.250;
    mt.mt[3] = 0.793;
    mt.mt[4] = 0.267;
    mt.mt[5] =-0.363;
    mt.sdr1[0] =335.0;
    mt.sdr1[1] = 46.0;
    mt.sdr1[2] =-113.0;
    mt.sdr2[0] = 186.0;
    mt.sdr2[1] = 49.0;
    mt.sdr2[2] =-68.0;
    mt.tax[0] = 2.32;
    mt.tax[1] = 2.0;
    mt.tax[2] = 261.0;
    mt.bax[0] = 0.38;
    mt.bax[1] = 16.0;
    mt.bax[2] = 351.0;
    mt.pax[0] =-2.70;
    mt.pax[1] = 74.0;
    mt.pax[2] =165.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[12], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposing 3\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------4------------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 4.9;
    mt.mt[0] =-2.270;
    mt.mt[1] = 0.058;
    mt.mt[2] = 2.210;
    mt.mt[3] = 0.079;
    mt.mt[4] =-0.737;
    mt.mt[5] = 0.246;
    mt.sdr1[0] =190.0;
    mt.sdr1[1] = 36.0;
    mt.sdr1[2] =-84.0;
    mt.sdr2[0] = 3.0;
    mt.sdr2[1] = 54.0;
    mt.sdr2[2] =-94.0;
    mt.tax[0] = 2.35;
    mt.tax[1] = 9.0;
    mt.tax[2] = 96.0;
    mt.bax[0] = 0.04;
    mt.bax[1] = 4.0;
    mt.bax[2] = 5.0;
    mt.pax[0] =-2.39;
    mt.pax[1] = 80.0;
    mt.pax[2] =253.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[18], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error dcomposing 4\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------5------------------------------------//
    mt.exp = 26.0 - 7.0;
    mt.mag = 6.5;
    mt.mt[0] = 0.745;
    mt.mt[1] =-0.036;
    mt.mt[2] =-0.709;
    mt.mt[3] =-0.242;
    mt.mt[4] =-0.048;
    mt.mt[5] = 0.208;
    mt.sdr1[0] =181.0;
    mt.sdr1[1] = 47.0;
    mt.sdr1[2] =114.0;
    mt.sdr2[0] = 328.0;
    mt.sdr2[1] = 48.0;
    mt.sdr2[2] = 67.0;
    mt.tax[0] = 0.82;
    mt.tax[1] = 73.0;
    mt.tax[2] = 166.0;
    mt.bax[0] =-0.05;
    mt.bax[1] = 17.0;
    mt.bax[2] = 344.0;
    mt.pax[0] =-0.77;
    mt.pax[1] = 1.0;
    mt.pax[2] = 74.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[24], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposing 5\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------6------------------------------------//
    mt.exp = 24.0 - 7.0;
    mt.mag = 5.3;
    mt.mt[0] = 0.700;
    mt.mt[1] =-0.883;
    mt.mt[2] = 0.183;
    mt.mt[3] = 0.260;
    mt.mt[4] = 0.289;
    mt.mt[5] = 0.712;
    mt.sdr1[0] =329.0;
    mt.sdr1[1] = 54.0;
    mt.sdr1[2] =141.0;
    mt.sdr2[0] = 85.0;
    mt.sdr2[1] = 59.0;
    mt.sdr2[2] = 43.0;
    mt.tax[0] = 1.01;
    mt.tax[1] = 51.0;
    mt.tax[2] = 300.0;
    mt.bax[0] = 0.24;
    mt.bax[1] = 39.0;
    mt.bax[2] = 113.0;
    mt.pax[0] =-1.25;
    mt.pax[1] = 3.0;
    mt.pax[2] = 206.;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[30], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposing 6\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------7------------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 5.0;
    mt.mt[0] = 3.150;
    mt.mt[1] =-2.470;
    mt.mt[2] =-0.676;
    mt.mt[3] = 1.650;
    mt.mt[4] = 1.880;
    mt.mt[5] =-1.470;
    mt.sdr1[0] =252.0;
    mt.sdr1[1] = 28.0;
    mt.sdr1[2] =111.0;
    mt.sdr2[0] = 49.0;
    mt.sdr2[1] = 64.0;
    mt.sdr2[2] = 79.0;
    mt.tax[0] = 4.08;
    mt.tax[1] = 69.0;
    mt.tax[2] = 297.0;
    mt.bax[0] = 0.01;
    mt.bax[1] = 10.0;
    mt.bax[2] = 54.0;
    mt.pax[0] =-4.08;
    mt.pax[1] = 18.0;
    mt.pax[2] = 147.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[36], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposing 7\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------8------------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 4.8;
    mt.mt[0] =-1.870;
    mt.mt[1] = 0.488;
    mt.mt[2] = 1.380;
    mt.mt[3] =-0.057;
    mt.mt[4] = 0.600;
    mt.mt[5] = 0.664;
    mt.sdr1[0] = 36.0;
    mt.sdr1[1] = 38.0;
    mt.sdr1[2] =-76.0;
    mt.sdr2[0] = 199.0;
    mt.sdr2[1] = 53.0;
    mt.sdr2[2] =-101.0;
    mt.tax[0] = 1.80;
    mt.tax[1] = 8.0;
    mt.tax[2] = 296.0;
    mt.bax[0] = 0.18;
    mt.bax[1] = 9.0;
    mt.bax[2] = 205.0;
    mt.pax[0] =-1.99;
    mt.pax[1] = 78.0;
    mt.pax[2] = 69.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[42], 1);
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decomposing 8\n");
        return EXIT_FAILURE;
    }
    //-----------------------------------9------------------------------------//
    mt.exp = 23.0 - 7.0;
    mt.mag = 5.0;
    mt.mt[0] =-3.540;
    mt.mt[1] = 1.040;
    mt.mt[2] = 2.510;
    mt.mt[3] =-0.474;
    mt.mt[4] = 1.640;
    mt.mt[5] = 1.530;
    mt.sdr1[0] = 47.0;
    mt.sdr1[1] = 38.0;
    mt.sdr1[2] =-63.0;
    mt.sdr2[0] = 195.0;
    mt.sdr2[1] = 56.0;
    mt.sdr2[2] =-109.0;
    mt.tax[0] = 3.66;
    mt.tax[1] = 10.0;
    mt.tax[2] = 299.0;
    mt.bax[0] = 0.45;
    mt.bax[1] = 16.0;
    mt.bax[2] =206.0;
    mt.pax[0] =-4.10;
    mt.pax[1] = 71.0;
    mt.pax[2] = 58.0;
    xscal = pow(10.0,mt.exp);
    cblas_dscal(n6, xscal, mt.mt, incx);
    cblas_dcopy(n6, mt.mt, 1, &mtTest[48], 1);
    // now copy this
    for (jmt=9; jmt<nmt; jmt++)
    {
        cblas_dcopy(n6, &mtTest[6*((jmt-9)%9)], 1, &mtTest[6*jmt], 1);
    }
    ierr = decompose(1, mt.mt);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error decompsoing 9\n");
        return EXIT_FAILURE;
    }
    //-----------------------------Test 'em all-------------------------------//
    ierr = decompose(nmt, mtTest);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "Error in bulk decomposition\n");
        return EXIT_FAILURE;
    }
    fprintf(stdout, "Success!\n");
    free(mtTest);
    return EXIT_SUCCESS;
}


