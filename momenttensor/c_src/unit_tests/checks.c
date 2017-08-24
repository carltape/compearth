#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#define CHKERR(t1, t2, tol, name, str)\
{\
   if (fabs(t1-t2) > tol){\
       printf("%s: %s %e %e\n", name, str, t1, t2);\
       return EXIT_FAILURE;\
   }\
}

int check_lam2nualpha(void);
int check_Udetcheck(void);
int check_lam2lune(void);
int check_CMTdecom(void);
int check_fangleSigned(void);
int check_CMT2faultpar(void);
int check_CMT2TT(void);
int check_normal2strdip(void);

int main(void)
{
    int ierr;

    ierr = check_lam2nualpha();
    if (ierr != 0){printf("failed lam2nualpha\n"); return EXIT_FAILURE;}
    printf("lam2nualpha was successful\n");

    ierr = check_fangleSigned();
    if (ierr != 0){printf("failed fangleSigned\n"); return EXIT_FAILURE;}
    printf("fangleSigned was successful\n");

    ierr = check_Udetcheck();
    if (ierr != 0){printf("failed Udetcheck\n"); return EXIT_FAILURE;}
    printf("udetcheck was successful\n");

    ierr = check_lam2lune(); 
    if (ierr != 0){printf("failed lam2lune\n"); return EXIT_FAILURE;}
    printf("lam2lune was succesful\n");

    ierr = check_CMTdecom();
    if (ierr != 0){printf("failed CMTdecom\n"); return EXIT_FAILURE;}
    printf("CMTdecom was successful\n");

    ierr = check_CMT2TT();
    if (ierr != 0){printf("failed CMT2TT\n"); return EXIT_FAILURE;}
    printf("CMT2TT was successful\n");

    ierr = check_CMT2faultpar();
    if (ierr != 0){printf("failed CMT2faultpar\n"); return EXIT_FAILURE;}
    printf("CMT2faultpar was successful\n");

    ierr = check_normal2strdip();
    if (ierr != 0){printf("failed noraml2strdip\n"); return EXIT_FAILURE;}
    printf("normal2strdip was successful\n");

    return EXIT_SUCCESS;
}
//============================================================================//
int check_fangleSigned(void)
{
    const char *fcnm = "check_fangleSigned\0"; 
    const double va[3] = {1, 0, 0};
    const double vb[3] = {0, 0, 3};
    const double vnor[3] = {0, -2, 0};
    int ierr;
    double stheta;
    stheta = compearth_eulerUtil_fangleSigned(3, va,vb,vnor, &ierr);
    if (fabs(stheta - 90.0) > 1.e-14 || ierr != 0)
    {
        printf("%s: failed test 1 %f %d\n", fcnm, stheta, ierr);
        return EXIT_FAILURE;
    }
    const double vnor2[3] =  {0, 2, 0};
    stheta = compearth_eulerUtil_fangleSigned(3, va,vb,vnor2, &ierr);
    if (fabs(stheta + 90.0) > 1.e-14 || ierr != 0)
    {
        printf("%s: failed test 2 %f %d\n", fcnm, stheta, ierr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int check_Udetcheck(void)
{
    const char *fcnm = "check_Udetcheck\0";
    const double U[9] = {0.6948,   0.3171,   0.9502,  // col 1
                         0.0344,   0.4387,   0.3816,  // col 2
                         0.7655,   0.7952,   0.1869}; // col 3
    const double Ur[9] = {0.6948,  0.3171,   0.9502,  // col 1
                         -0.0344, -0.4387,  -0.3816,  // col 2 negated
                          0.7655,  0.7952,   0.1869};
    double Un[9], Un2[9];
    int i, ierr;
    ierr = compearth_Udetcheck(1, U, Un);
    if (ierr != 0)
    {
        printf("%s: Error in Udetcheck1\n", fcnm);
        return EXIT_FAILURE;
    }
    ierr = compearth_Udetcheck(1, Un, Un2);
    if (ierr != 0)
    {   
        printf("%s: Error in Udetcheck2\n", fcnm);
        return EXIT_FAILURE;
    }
    for (i=0; i<9; i++)
    {
        if (fabs(Ur[i] - Un[i]) > 1.e-12)
        {
            printf("%s: Udetcheck failed; %f %f\n", fcnm, Ur[i], Un[i]);
            return EXIT_FAILURE;
        }
        if (fabs(Un[i] - Un2[i]) > 1.e-12)
        {
            printf("%s: Udetcheck2 failed; %f %f\n", fcnm, Un[i], Un2[i]);
            return EXIT_FAILURE;
        }
    } 
    return EXIT_SUCCESS;
}
//============================================================================//
int check_lam2nualpha(void)
{ 
    const char *fcnm = "check_lam2nualpha\0";
    // example from TT2013, App A.
    double lam[3] = { 8.802, 2.584, -1.851};
    double lamcheck[3], nu, alpha, mag; 
    const double nu1 = 0.371745072651417;
    const double alpha1 = 80.365019327257144;
    int ierr;
    ierr = compearth_lam2nualpha(1, lam, &nu, &alpha);
    if (ierr != 0)
    {
        printf("%s: error calling lam2nualpha\n", fcnm);
        return EXIT_FAILURE;
    }
    CHKERR(nu, nu1, 1.e-10, fcnm, "error computing nu 1");
    CHKERR(alpha, alpha1, 1.e-10, fcnm, "error computing alpha 1"); 
    ierr = compearth_nualpha2lam(1, &nu, &alpha, lamcheck);
    if (ierr != 0)
    {
        printf("%s: error caling nualpha2lam\n", fcnm);
        return EXIT_FAILURE;
    } 
    mag = cblas_dnrm2(3, lam, 1);
    lam[0] = lam[0]/mag;
    lam[1] = lam[1]/mag;
    lam[2] = lam[2]/mag;
    CHKERR(lam[0], lamcheck[0], 1.e-10, fcnm, "error computing lam 1");
    CHKERR(lam[1], lamcheck[1], 1.e-10, fcnm, "error computing lam 2");
    CHKERR(lam[2], lamcheck[2], 1.e-10, fcnm, "error computing lam 3");
    //lam / norm(lam)
    return EXIT_SUCCESS;
}
//============================================================================//
int check_CMT2faultpar(void)
{
    const double M[6] = {3.108304932835845, 3.044425632830430, 3.382269434333724,
                     -4.855033301709626,-1.949280336439431, 1.110527600460120};
    double nu, alpha, N1[3], N2[3], lam[3];
    const double nuRef = 0.371745072651417;
    const double alphaRef = 80.365019327257144;
    const double N1ref[3] = {-0.050351975426048, 0.930119048707825, 0.363790095799137};
    const double N2ref[3] = {-0.977452142939121, 0.046467470073259, 0.205980781843140};
    const double lamref[3] = { 8.802000000000001, 2.584000000000001, -1.851000000000002};
    int ierr;
    int nmt = 1;
    ierr = compearth_CMT2faultpar(nmt, M, 
                                  &nu, &alpha, N1, N2, lam);
    if (ierr != 0)
    {
        fprintf(stderr, "Error calling CMT2faultpar\n");
        return EXIT_FAILURE;
    }
    CHKERR(nu,    nuRef,    1.e-10, __func__, "error computing nu");
    CHKERR(alpha, alphaRef, 1.e-10, __func__, "error computing alpha");
    CHKERR(N1[0], N1ref[0], 1.e-10, __func__, "error checking N1[0]");
    CHKERR(N1[1], N1ref[1], 1.e-10, __func__, "error checking N1[1]");
    CHKERR(N1[2], N1ref[2], 1.e-10, __func__, "error checking N1[2]");
    CHKERR(N2[0], N2ref[0], 1.e-10, __func__, "error checking N2[0]");
    CHKERR(N2[1], N2ref[1], 1.e-10, __func__, "error checking N2[1]");
    CHKERR(N2[2], N2ref[2], 1.e-10, __func__, "error checking N2[2]");
    CHKERR(lam[0], lamref[0], 1.e-10, __func__, "error checking lam[0]");
    CHKERR(lam[1], lamref[1], 1.e-10, __func__, "error checking lam[1]");
    CHKERR(lam[2], lamref[2], 1.e-10, __func__, "error checking lam[2]");
    return EXIT_SUCCESS;
}
//============================================================================//
int check_CMT2TT(void)
{
    const char *fcnm = "check_CMT2TT\0";
    const double M[6] = {3.108304932835845, 3.044425632830430,
                         3.382269434333724,-4.855033301709626,
                        -1.949280336439431, 1.110527600460120};
    double U[9], lam[3];
    bool ldisplay;
    double gamma, delta, M0, kappa, sigma, theta, thetadc;
    double K[3], N[3], S[3];
    const double gamma1 = -5.519441015228609;
    const double delta1 = 36.032871051674995;
    const double M01 = 6.617343160211657;
    const double kappa1 = 69.492561954605137;
    const double theta1 = 88.144957350922951;
    const double sigma1 = -79.788954797632840;
    const double thetadc1 = 36.396475670817210;
    const double KR1[3] = {-0.350328975576259, 0.936626718000127, 0};
    const double NR1[3] = {0.936135854045508, 0.350145376429380, 0.032370945856051};
    const double SR1[3] = {-0.032265105849630, 0.177200864349715, -0.983645676358224};
    const double lam1[3] = {8.802000000000000, 2.584000000000003, -1.851000000000001};
    const double UR1[9] = {-0.639133135365464,-0.372890102888132, 0.672652812709492,
                            0.350155145207092,-0.919781533321278,-0.177181560118875,
                            0.684762885649414, 0.122290237260530, 0.718432243365964};
    // Test 2; problem child of horizontal fault
    const double M2[6] = {0, 0, 0, -sqrt(3)/2., .5, 0};
    const double gamma2 = 0.0; const double delta2 = 0.0;
    const double M02 = 1.0; const double kappa2 = 300.0;
    const double theta2 = 90.0; const double sigma2 = 90.0;
    const double thetadc2 = 1.207418269725733e-06;
    const double KR2[3] = {-0.500000000000000, -0.866025403784439, 0};
    const double NR2[3] = {-0.866025403784439, 0.500000000000000,  0};
    const double SR2[3] = {0, 0, 1};
    const double lam2[3] = {1, 0, -1};
    const double UR2[9] = {-0.612372435695794,0.353553390593274,0.707106781186547,
                            0.500000000000000,0.866025403784439,0,
                           -0.612372435695794,0.353553390593274,-0.707106781186547};
    
    int i, ierr;
    ldisplay = false;
    ierr = compearth_CMT2TT(1, M, ldisplay,
                            &gamma, &delta, &M0, &kappa, &theta, &sigma,
                            K, N, S, &thetadc,
                            lam, U);
    if (ierr != 0)
    {
        printf("%s: Error calling CMT2TT test 1\n", fcnm);
        return EXIT_FAILURE;
    }
    CHKERR(gamma, gamma1, 1.e-10, fcnm, "error computing gamma 1"); 
    CHKERR(delta, delta1, 1.e-10, fcnm, "error computing delta 1");
    CHKERR(M0, M01, 1.e-10, fcnm, "error computing M0 1");
    CHKERR(kappa, kappa1, 1.e-10, fcnm, "error computing kappa 1");
    CHKERR(theta, theta1, 1.e-10, fcnm, "error computing theta 1");
    CHKERR(sigma, sigma1, 1.e-10, fcnm, "error computing sigma 1");
    CHKERR(thetadc, thetadc1, 1.e-10, fcnm, "error computing thetadc 1");
    for (i=0; i<3; i++)
    {
        if (fabs(KR1[i] - K[i]) > 1.e-10)
        {
            printf("%s: Error computing K 1\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(NR1[i] - N[i]) > 1.e-10)
        {
            printf("%s: Error computing N 1\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(SR1[i] - S[i]) > 1.e-10)
        {
            printf("%s: Error computing S 1\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(lam1[i] - lam[i]) > 1.e-10)
        {
            printf("%s: Error computing lam 1\n", fcnm);
            return EXIT_FAILURE;
        }
    }
    for (i=0; i<9; i++)
    {
        if (fabs(U[i] - UR1[i]) > 1.e-10)
        {
            printf("%s: Error computing U 1\n", fcnm);
            return EXIT_FAILURE;
        }
    }

    ierr = compearth_CMT2TT(1, M2, ldisplay,
                            &gamma, &delta, &M0, &kappa, &theta, &sigma,
                            K, N, S, &thetadc,
                            lam, U);
    if (ierr != 0)
    {
        printf("%s: Error calling CMT2TT horizontal fault\n", fcnm);
        return EXIT_FAILURE;
    }
    CHKERR(gamma, gamma2, 1.e-10, fcnm, "error computing gamma 2"); 
    CHKERR(delta, delta2, 1.e-10, fcnm, "error computing delta 2");
    CHKERR(M0, M02, 1.e-10, fcnm, "error computing M0 2");
    CHKERR(kappa, kappa2, 1.e-10, fcnm, "error computing kappa 2");
    CHKERR(theta, theta2, 1.e-10, fcnm, "error computing theta 2");
    CHKERR(sigma, sigma2, 1.e-10, fcnm, "error computing sigma 2");
    CHKERR(thetadc, thetadc2, 1.e-10, fcnm, "error computing thetadc 2");
    for (i=0; i<3; i++)
    {   
        if (fabs(KR2[i] - K[i]) > 1.e-10)
        {
            printf("%s: Error computing K 2\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(NR2[i] - N[i]) > 1.e-10)
        {
            printf("%s: Error computing N 2\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(SR2[i] - S[i]) > 1.e-10)
        {
            printf("%s: Error computing S 2\n", fcnm);
            return EXIT_FAILURE;
        }
        if (fabs(lam2[i] - lam[i]) > 1.e-10)
        {
            printf("%s: Error computing lam 2\n", fcnm);
            return EXIT_FAILURE;
        }
    }
    for (i=0; i<9; i++)
    {
        if (fabs(U[i] - UR2[i]) > 1.e-10)
        {
            printf("%s: Error computing U 2\n", fcnm);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int check_CMTdecom(void)
{
    const char *fcnm = "check_CMTdecom\0";
    const double M[6] = {3.108304932835845, 3.044425632830430,
                         3.382269434333724,-4.855033301709626,
                        -1.949280336439431, 1.110527600460120};
    double U[9], lam[3];
    double lamRef1[3] = {8.8020000000e+00, 2.5840000000e+00, -1.8510000000e+00};
    double lamRef2[3] = {-1.851000000000002, 2.584000000000001, 8.802000000000001};
    double *lamRef3 = lamRef1;
    double *lamRef4 = lamRef2;
    double URef1[9] = {-0.672652812709492, 0.639133135365464, 0.372890102888132,
                        0.177181560118875,-0.350155145207092, 0.919781533321279,
                        0.718432243365964, 0.684762885649414, 0.122290237260530};
    double URef2[9] = { 0.718432243365964, 0.684762885649414, 0.122290237260530,
                       -0.177181560118875, 0.350155145207092,-0.919781533321279,
                       -0.672652812709492, 0.639133135365464, 0.372890102888132};
    double *URef3 = URef1;
    double *URef4 = URef2;
    int i, ierr;
    ierr = compearth_CMTdecom(1,  M, 1, lam, U);
    if (ierr != 0)
    {
        printf("%s: Error 1 calling CMTdecom\n", fcnm);
        return EXIT_FAILURE;
    } 
    for (i=0; i<3; i++)
    {
        if (fabs(lam[i] - lamRef1[i]) > 1.e-13)
        {
            printf("%s: lam 1 is wrong %e %e\n", fcnm, lam[i], lamRef1[i]);
            return EXIT_FAILURE;
        }
    }
    for (i=0; i<9; i++)
    {
        if (fabs(U[i] - URef1[i]) > 1.e-13)
        {
            printf("%s: U 1 is wrong %d %e %e\n", fcnm, i, U[i], URef1[i]);
            return EXIT_FAILURE;
        }
    }

    ierr = compearth_CMTdecom(1,  M, 2, lam, U); 
    if (ierr != 0)
    {   
        printf("%s: Error 2 calling CMTdecom\n", fcnm);
        return EXIT_FAILURE;
    }   
    for (i=0; i<3; i++)
    {
        if (fabs(lam[i] - lamRef2[i]) > 1.e-13)
        {
            printf("%s: lam 2 is wrong %e %e\n", fcnm, lam[i], lamRef2[i]);
            return EXIT_FAILURE;
        }
    }
    for (i=0; i<9; i++)
    {
        if (fabs(U[i] - URef2[i]) > 1.e-13)
        {
            printf("%s: U 2 is wrong %d %e %e\n", fcnm, i, U[i], URef2[i]);
            return EXIT_FAILURE;
        }
    }

    ierr = compearth_CMTdecom(1,  M, 3, lam, U); 
    if (ierr != 0)
    {   
        printf("%s: Error calling CMTdecom 3\n", fcnm);
        return EXIT_FAILURE;
    }   
    for (i=0; i<3; i++)
    {   
        if (fabs(lam[i] - lamRef3[i]) > 1.e-13)
        {
            printf("%s: lam 3 is wrong %e %e\n", fcnm, lam[i], lamRef3[i]);
            return EXIT_FAILURE;
        }
    }   
    for (i=0; i<9; i++)
    {
        if (fabs(U[i] - URef3[i]) > 1.e-13)
        {
            printf("%s: U 3 is wrong %d %e %e\n", fcnm, i, U[i], URef3[i]);
            return EXIT_FAILURE;
        }
    }

    ierr = compearth_CMTdecom(1,  M, 4, lam, U); 
    if (ierr != 0)
    {   
        printf("%s: Error calling CMTdecom 4\n", fcnm);
        return EXIT_FAILURE;
    }   
    for (i=0; i<3; i++)
    {   
        if (fabs(lam[i] - lamRef4[i]) > 1.e-13)
        {
            printf("%s: lam 4 is wrong %e %e\n", fcnm, lam[i], lamRef4[i]);
            return EXIT_FAILURE;
        }
    }   
    for (i=0; i<9; i++)
    {   
        if (fabs(U[i] - URef4[i]) > 1.e-13)
        {
            printf("%s: U 4 is wrong %d %e %e\n", fcnm, i, U[i], URef4[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int check_lam2lune(void)
{
    const char *fcnm = "check_lam2lune\0";
    const double lamRef1[6] = {
        -7.9621992143e+15,-8.2644837610e+15,-8.2644837610e+15,
        -7.5915701171e+15,-8.4370629711e+15,-8.4370629711e+15};
    double *gamma0, *gamma, *delta, *delta0, *lam, *M0, *M00;
    double *thetadc, *lamdev, *lamiso;
    int i, ib, ierr, ig, indx;
    const int ng = 100;
    const int nb = 100;
    int nmt = ng*nb;
    const double dg = (30.0 - -30.0)/(double) (ng - 1);
    const double db = (89.0 - -89.0)/(double) (nb - 1);
    gamma0 = (double *) calloc((size_t) nmt, sizeof(double));
    delta0 = (double *) calloc((size_t) nmt, sizeof(double));
    M00    = (double *) calloc((size_t) nmt, sizeof(double));
    gamma  = (double *) calloc((size_t) nmt, sizeof(double));
    delta  = (double *) calloc((size_t) nmt, sizeof(double));
    M0     = (double *) calloc((size_t) nmt, sizeof(double)); 
    lam  = (double *) calloc((size_t) (3*nmt), sizeof(double)); 
    for (ig=0; ig<ng; ig++)
    {
        for (ib=0; ib<nb; ib++)
        {
            indx = ig*nb + ib;
            gamma0[indx] =-30.0 + (double) ig*dg;
            delta0[indx] =-89.0 + (double) ib*db; // use latitude instead of colat
            M00[indx] = 1.e16;
        }
    }
    compearth_lune2lam(nmt, gamma0, delta0, M00, lam);
    // Verify the first two computations are right
    for (i=0; i<6; i++)
    {
        if (fabs(lam[i] - lamRef1[i])/1.e16 > 1.e-10)
        {
            printf("%s: error in lune2lam: %.10e %.10e\n",
                   fcnm, lam[i], lamRef1[i]);
            return EXIT_FAILURE;
        }
    }
    // Now perform the inverse operation
    thetadc = NULL;
    lamdev = NULL;
    lamiso = NULL;
    ierr = compearth_lam2lune(nmt, lam, 
                              gamma, delta, M0,
                              thetadc, lamdev, lamiso);
    if (ierr != 0)
    {
        printf("%s: error calling lam2lune\n", fcnm);
        return EXIT_FAILURE;
    }
    // See how close i was
    for (i=0; i<nmt; i++)
    {
        if (fabs(gamma[i] - gamma0[i]) > 1.e-12)
        {
            printf("%s: gamma different %e %e\n", fcnm, gamma[i], gamma0[i]);
        }
        if (fabs(delta[i] - delta0[i]) > 1.e-12)
        {
            printf("%s: delta different %e %e\n", fcnm, delta[i], delta[i]);
        }
        if (fabs(M0[i] - M00[i])/1.e15 > 1.e-12)
        {
           printf("%s: M0 different %e %e\n", fcnm, M0[i], M00[i]);
        }
    }
    free(gamma0);
    free(delta0);
    free(M00);
    free(gamma);
    free(delta);
    free(M0);
    free(lam);
    return EXIT_SUCCESS;
}

int check_normal2strdip(void)
{
    const double Xin[3] = {-0.171010071662835, 0.969846310392954,
                            0.173648177666930};
    double Xout[2];
    compearth_normal2strdip(1, Xin, Xout);
    CHKERR(Xout[0], 350.0, 1.e-10, __func__, "error computing Xout[0]"); 
    CHKERR(Xout[1], 80.0,  1.e-10, __func__, "error computing Xout[1]");
    return EXIT_SUCCESS;
}
