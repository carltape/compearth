#ifndef COMPEARTH_H__
#define COMPEARTH_H__ 1
#include <stdbool.h>
#include <math.h>
#include "compearth_config.h"
#include "compearth_enum.h"
#include "compearth_constants.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Angle between moment tensors */
int compearth_angleMT(const int n,
                      const double *__restrict__ M1in,
                      const double *__restrict__ M2in,
                      double *__restrict__ theta);
/* Auxiliary plane */
int compearth_auxiliaryPlane(const int nmt,
                             const double *__restrict__ s1, 
                             const double *__restrict__ d1, 
                             const double *__restrict__ r1, 
                             double *__restrict__ s2, 
                             double *__restrict__ d2, 
                             double *__restrict__ r2);
/* Converts u to lune colatitude beta */
int compearth_u2beta(const int n, const double *__restrict__ u,
                     double *__restrict__ beta);
int compearth_u2beta_optimize(const int n,
                              const int maxit,
                              const int linvType,
                              const double *__restrict__ u,
                              const double tol,
                              double *__restrict__ beta);
/* Converts lune colatitude beta to u */
void compearth_beta2u(const int n, const double *__restrict__ beta,
                      double *__restrict__ u);
/* Decompose a CMT */
int compearth_CMTdecom(const int nmt, const double *__restrict__ M,
                       const int isort,
                       double *__restrict__ lam,
                       double *__restrict__ U);
/* Decompose moment tensor into deviatori and isotropic parts */
int compearth_CMTdecomIso(const int n, const double *__restrict__ M,
                          double *__restrict__ Miso,
                          double *__restrict__ Mdev,
                          double *__restrict__ trM);
/* Get unit normals to fault planes from moment tensor M */
int compearth_CMT2faultpar(const int nmt,
                          const double *__restrict__ M,
                          double *__restrict__ nu,
                          double *__restrict__ alpha,
                          double *__restrict__ N1,
                          double *__restrict__ N2,
                          double *__restrict__ lam);
/* Compute M0 from CMT */
int compearth_CMT2m0(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ M0);
/* Compute Mw from CMT */
int compearth_CMT2mw(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ Mw);
/* Compute angle between moment tensors */
int compearth_CMT2omega(const int nmt1, const double *__restrict__ M1, 
                        const int nmt2, const double *__restrict__ M2, 
                        double *__restrict__ omega);
/* CMT to Tape and Tape */
int compearth_CMT2TT(const int nmt, const double *__restrict__ Min,
                     const bool ldisplay,
                     double *__restrict__ gamma,
                     double *__restrict__ delta,
                     double *__restrict__ M0,
                     double *__restrict__ kappa,
                     double *__restrict__ theta,
                     double *__restrict__ sigma,
                     double *__restrict__ K, double *__restrict__ N,
                     double *__restrict__ S, double *__restrict__ thetadc,
                     double *__restrict__ lam, double *__restrict__ U);
/* Converts between moment tensors in different coordinate systems */
int compearth_convertMT(const int nmt,
                        const enum compearthCoordSystem_enum i1in,
                        const enum compearthCoordSystem_enum i2in,
                        const double *__restrict__ M,
                        double *__restrict__ Mout);
/* Computes the signed angle between two vectors */
double compearth_eulerUtil_fangleSigned(const int n,
                                        const double *__restrict__ va, 
                                        const double *__restrict__ vb, 
                                        const double *__restrict__ vnor,
                                        int *ierr);
/* Computes rotation matrices about axis for given angles */
int compearth_eulerUtil_rotmat(const int n, const double *__restrict__ xdeg,
                               const int ixyz, double *__restrict__ R);
/* Compute rotation matrices about axis for given angles */
int compearth_eulerUtil_rotmatGen(const int n,
                                  const double *__restrict__ v,
                                  const double *__restrict__ xi,
                                  double *__restrict__ U);
/* Converts lune longitude gamma to v */
void compearth_gamma2v(const int n, const double *__restrict__ gamma,
                       double *__restrict__ v);
/* Converts eigenvalues to nu and alpha */
int compearth_lam2nualpha(const int n, const double *__restrict__ lam,
                          double *__restrict__ nu, double *__restrict__ alpha);
/* Converts eigenvalues to phi and zeta */
int compearth_lam2phizeta(const int n, const double *__restrict__ lam,
                          double *__restrict__ phi,
                          double *__restrict__ zeta);
/* Converts lambda (eigenvalues) to lune values */
int compearth_lam2lune(const int nmt, const double *__restrict__ lam,
                       double *__restrict__ gamma,
                       double *__restrict__ delta,
                       double *__restrict__ M0, 
                       double *__restrict__ thetadc,
                       double *__restrict__ lamdev,
                       double *__restrict__ lamiso);
/* Sorts eigenvalues */
int compearth_lamsort(const int n, const double *__restrict__ lam,
                      double *__restrict__ lamSort);
/* Convert lune coordinates to eigenvalues */
void compearth_lune2lam(const int n,
                        const double *__restrict__ gamma,
                        const double *__restrict__ delta,
                        const double *__restrict__ M0in,
                        double *__restrict__ lam);
/* Convert lune to rect coordinates */
void compearth_lune2rect(const int ng, double *__restrict__ gamma,
                         const int nd, double *__restrict__ delta,
                         double *__restrict__ v,
                         double *__restrict__ w);
/* cartesian to spherical coordinate conversion */
void compearth_matlab_cart2sph(const int n,
                               const double *__restrict__ x,
                               const double *__restrict__ y,
                               const double *__restrict__ z,
                               double *__restrict__ theta,
                               double *__restrict__ phi,
                               double *__restrict__ r);
/* angle between two vectors */
double compearth_matlab_fangle(const int n,
                               const double *__restrict__ va, 
                               const double *__restrict__ vb);
/* M0 to half duration */
int compearth_m02hdur(const int nm, const double *__restrict__ M0,
                      double *__restrict__ hdur);
/* Compute Mw from M0 */
int compearth_m02mw(const int nm, const enum magType_enum imag,
                    const double *__restrict__ M0, 
                    double *__restrict__ Mw);
/* Compute M0 from Mw */
int compearth_mw2m0(const int nm, const enum magType_enum imag,
                    const double *__restrict__ Mw, 
                    double *__restrict__ M0);
/* Convert vector moment tensors to matrix moment tensors */
void compearth_Mvec2Mmat(const int n,
                         const double *__restrict__ Min,
                         const int itype,
                         double *__restrict__ Mout);
/* Converts fault normal vectors to strike and dip angles */
void compearth_normal2strdip(const int n,
                             const double *__restrict__ Xin,
                             double *__restrict__ Xout);
/* Norm of matrix stored as an array */
int compearth_normMat(const int n,
                      const double *M, 
                      const enum ceNormType_enum Lnorm,
                      const double p,
                      double *mnorm);
/* Norm of a moment tensor */
int compearth_normMT(const int n,
                     const double *M, 
                     const enum ceNormType_enum Lnorm,
                     const double p,
                     double *mnorm);
/* Converts (nu,alpha) to unit lambda vector */
int compearth_nualpha2lam(const int n,
                          const double *__restrict__ nu,
                          const double *__restrict__ alpha,
                          double *__restrict__ lam);
/* Converts (v,w) rect coordinates to lune coordinates */
int compearth_rect2lune(const int nv, const double *__restrict__ v,
                        const int nw, const double *__restrict__ w,
                        double *__restrict__ gamma,
                        double *__restrict__ delta);
/* Lune coordinates to eigenvalue triple */
void compearth_tape2015Eqn7(const double beta, const double gamma,
                            double lam[3]);
/* Convert a tape and tape parameterized moment tensor to a moment tensor */
int compearth_TT2CMT(const int nmt,
                     const double *__restrict__ gamma,
                     const double *__restrict__ delta,
                     const double *__restrict__ M0,
                     const double *__restrict__ kappa,
                     const double *__restrict__ theta,
                     const double *__restrict__ sigmaIn,
                     double *__restrict__ M,
                     double *__restrict__ lam,
                     double *__restrict__ U);
/*
int compearth_tt2cmt(const double gamma,
                     const double delta,
                     const double M0,
                     const double kappa,
                     const double theta,
                     const double sigmaIn,
                     double M[6], double lam[3], double U[9]);
*/
/* Converts theta dip angle to h */
void compearth_theta2h(const int n, const double *__restrict__ theta,
                        double *__restrict__ h);
/* U to plunge and azimuth */
int compearth_U2pa(const int nmt, const double *__restrict__ U,
                   double *__restrict__ pl1,
                   double *__restrict__ az1,
                   double *__restrict__ pl2,
                   double *__restrict__ az2,
                   double *__restrict__ pl3,
                   double *__restrict__ az3);
/* Ensures det >= 0 for rotation matrices. */
int compearth_Udetcheck(const int n,
                        const double *__restrict__ Uin,
                        double *__restrict__ Uout);
/* Orthogonalizes a 3 x 3 matrix. */
int compearth_Uorth(const int n, const enum ceOrthoType_enum itype,
                    const double *__restrict__ Uin,
                    double *__restrict__ Uout,
                    double *__restrict__ dtUin,
                    double *__restrict__ dtUout);
/* Transform moment tensor with transformation matrix T */
int compearth_transform_mat(const int nmt,
                            const double *__restrict__ T,
                            const double *__restrict__ Min,
                            double *__restrict__ Mout);
/* Converts gamma to lune longitude v */
void compearth_v2gamma(const int n, const double *__restrict__ v,
                       double *__restrict__ gamma);
/* Maps theta dip angle to h */
void compearth_h2theta(const int n, const double *__restrict__ h,
                       double *__restrict__ theta);
/* xyz to latitude, longtiude, radius */
void compearth_xyz2tp(const int n,
                      const double *__restrict__ x,
                      const double *__restrict__ y,
                      const double *__restrict__ z,
                      double *__restrict__ ph, 
                      double *__restrict__ th, 
                      double *__restrict__ rho);
/* xyz to lat lon */
void compearth_xyz2latlon(const int n,
                          const double *__restrict__ x,
                          const double *__restrict__ y,
                          const double *__restrict__ z,
                          double *__restrict__ lat,
                          double *__restrict__ lon); 


/*! Standard decomposition */
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

#ifdef COMPEARTH_PRIVATE_DET3X3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern double det3x3ColumnMajor(const double *__restrict__ A);
#endif
#ifdef COMPEARTH_PRIVATE_GEMV3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern void gemv3_colMajorNoTrans(const double *__restrict__ A,
                                  const double *__restrict__ x,
                                  double *__restrict__ y);
#endif
#ifdef COMPEARTH_PRIVATE_GEM3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern void gemm3_colMajorNoTransNoTrans(const double *__restrict__ A,
                                         const double *__restrict__ B,
                                         double *__restrict__ C);
#endif
#ifdef COMPEARTH_PRIVATE_GEMT3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern void gemm3_colMajorNoTransTrans(const double *__restrict__ A,
                                       const double *__restrict__ B,
                                       double *__restrict__ C);
#endif
#ifdef COMPEARTH_PRIVATE_CROSS3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern void cross3(const double *__restrict__ a, 
                   const double *__restrict__ b,
                   double *__restrict__ c);
#endif
#ifdef COMPEARTH_PRIVATE_DOT3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern double dot3(const double *__restrict__ a, const double *__restrict__ b);
#endif
#ifdef COMPEARTH_PRIVATE_NORM3
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern double norm3(const double *__restrict__ a);
#endif
#ifdef COMPEARTH_PRIVATE_WRAP360
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern double wrap360(const double lon);
#endif


#ifdef COMPEARTH_PRIVATE_MOD
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern double mod(const double x, const double y);
#endif

#ifdef COMPEARTH_PRIVATE_ANTIPODE
#ifdef _OPENMP
#pragma omp declare simd
#endif
extern void antipode(const double latIn, const double lonIn,
                     double *latOut, double *lonOut, const bool isDeg);
#endif

#ifdef COMPEARTH_PRIVATE_UPDOWN_ARGSORT3
extern int argsort3_upDown(const double *__restrict__ x,
                           const bool lascend,
                           int *__restrict__ perm);
#endif
#ifdef COMPEARTH_PRIVATE_UPDOWN_ABS_ARGSORT3
extern int argsort3_absUpDown(const double *__restrict__ x,
                              const bool lascend,
                              int *__restrict__ perm);
#endif
#ifdef COMPEARTH_PRIVATE_ARGSORT3
extern int argsort3(const double *__restrict__ x, int *__restrict__ perm);
#endif

#ifdef __cplusplus
}
#endif

#endif
