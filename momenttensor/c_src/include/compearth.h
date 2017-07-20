#ifndef COMPEARTH_H__
#define COMPEARTH_H__ 1
#include <stdbool.h>
#include <math.h>
#include "compearth_config.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

enum compearthCoordSystem_enum
{
    USE = 1,    /*!< Up, south, east (G-CMT) */
    NED = 2,    /*!< North, east, down (Aki and Richards, 1980 pg 118) */
    NWU = 3,    /*!< North, west, up */
    ENU = 4,    /*!< East, north, up */
    SEU = 5     /*!< South, east, up */
};

enum magType_enum
{
    KANAMORI_1978 = 1, /*!< Mw frrom Kanamori Mw = (2/3)*log10(M0) + k; */
    HARVARD_CMT = 2    /*!< Mw from Harvard CMT: (2/3)*(log10(M0) - 16.1) */
};

#ifdef COMPEARTH_USE_ISCL
#include "iscl/array/array.h"
#else
enum normType_enum
{
    TWO_NORM = 2,              /*!< \$ L_2 = \sqrt{\sum_i x_i^2} \$ norm */
    ONE_NORM = 1,              /*!< \$ L_1 = \sum_i |x_i| \$ norm */
    P_NORM = 3,                /*!< \$ L_p 
                                    = \left (
                                        \sum_i |x_i|^p \right )^{1/p}
                                      \right ) \$ norm */
    INFINITY_NORM = 4,          /*!< \$ L_\infty = max |x| \$ */
    NEGATIVE_INFINITY_NORM = 5  /*!< \$ L_{-\infty} = min |x| \$ */
};
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/* Angle between moment tensors */
int compearth_angleMT(const int n,
                      const double *__restrict__ M1in,
                      const double *__restrict__ M2in,
                      double *__restrict__ theta);
/* Converts u to lune colatitude beta */
int compearth_u2beta(const int n,
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
/* Compute M0 from CMT */
int compearth_CMT2m0(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ M0);
/* Compute Mw from CMT */
int compearth_CMT2mw(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ Mw);
/* Converts between moment tensors in different coordinate systems */
int compearth_convertMT(const enum compearthCoordSystem_enum i1in,
                        const enum compearthCoordSystem_enum i2in,
                        const double *__restrict__ M,
                        double *__restrict__ Mout);
/* Computes rotation matrices about axis for given angles */
int compearth_eulerUtil_rotmat(const int n, const double *__restrict__ xdeg,
                               const int ixyz, double *__restrict__ R);
/* Compute rotation matrices about axis for given angles */
int compearth_eulerUtil_rotmatGen(const int n,
                                  const double *__restrict__ v,
                                  const double *xi,
                                  double *__restrict__ U);
/* Converts lune longitude gamma to v */
void compearth_gamma2v(const int n, const double *__restrict__ gamma,
                       double *__restrict__ v);
/* Converts eigenvalues to phi and zeta */
int compearth_lam2phizeta(const int n, const double *__restrict__ lam,
                          double *__restrict__ phi,
                          double *__restrict__ zeta);
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
                      const enum normType_enum Lnorm,
                      const double p,
                      double *mnorm);
/* Norm of a moment tensor */
int compearth_normMT(const int n,
                     const double *M, 
                     const enum normType_enum Lnorm,
                     const double p,
                     double *mnorm);
/* Converts (v,w) rect coordinates to lune coordinates */
int compearth_rect2lune(const int nv, const double *__restrict__ v,
                        const int nw, const double *__restrict__ w,
                        double *__restrict__ gamma,
                        double *__restrict__ delta);
/* Lune coordinates to eigenvalue triple */
void compearth_tape2015Eqn7(const double beta, const double gamma,
                            double lam[3]);
/* Convert a tape and tape parameterized moment tensor to a moment tensor */
int compearth_tt2cmt(const double gamma,
                     const double delta,
                     const double M0,
                     const double kappa,
                     const double theta,
                     const double sigmaIn,
                     double M[6], double lam[3], double U[9]);
/* Converts theta dip angle to h */
void compearth_theta2h(const int n, const double *__restrict__ theta,
                        double *__restrict__ h);
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
#ifdef __cplusplus
}
#endif

#endif
