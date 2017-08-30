#ifndef _cmopad_h__
#define _cmopad_h__ 1
#include <stdbool.h>

struct cmopad_struct
{
    double M[3][3];               /*!< Full moment tensor */
    double M_devi[3][3];          /*!< Deviatoric moment tensor (i.e. no
                                       isotropic component) */
    double M_iso[3][3];           /*!< Isotropic moment tensor */
    double DC[3][3];              /*!< Best fitting double couple moment
                                       tensor */
    double CLVD[3][3];            /*!< Compensated linear vector dipole moment
                                       tensor */
    double eigenvectors[3][3];    /*!< Internal workspace */
    double rotation_matrix[3][3]; /*!< Internal workspace */
    double eigenvalues[3];        /*!< Workspace holds eigenvalues */
    double eig_pnt[3];            /*!< Eigenvalues for pressure, null,
                                       and tension */
    double null_axis[3];          /*!< Null axis eigenvector */
    double t_axis[3];             /*!< Tension axis eigenvector */
    double p_axis[3];             /*!< Pressure axis eigenvector */
    double fp1[3];                /*!< Strike, dip, and rake of fault
                                       describing plane 1 (degrees) */
    double fp2[3];                /*!< Strike, dip, and rake of fault
                                       describing plane 2 (degrees) */
    double seismic_moment;        /*!< Scalar moment (Nm) */
    double moment_magnitude;      /*!< Moment magnitude (e.g. Kanamori 1977) */

    double DC_percentage;         /*!< Double couple percentage [0,100] */
    double ISO_percentage;        /*!< Isotrpic percentage [0,100] */
    double DEV_percentage;        /*!< Percent deviatoric [0,100] 
                                       This is 100 - DC_percentage */
    double CLVD_percentage;       /*!< Percentage of DC part of the moment tensor
                                       tensor.  This is:
                                       100 - ISO_percentage - DC_percentage */

    int plot_clr_order;           /*!< Pertinent to ordering the nodal plane
                                        ordering with regards to plotting */
    char pad[4];
};

enum cmopad_basis_enum
{
    NED = 1,   /*!< North, East, Down - like in Jost and Herrmann */
    USE = 2,   /*!< Up, South, East - like in Global CMT */
    XYZ = 3,   /*!< North, Esat, Up */
    NWU = 4    /*!< North, West, Up */
};

#ifdef __cplusplus
extern "C" 
{
#endif

void cmopad_findStrikeDipRake(double rot_mat[3][3],
                              double *alpha, double *beta, double *gamma);
void cmopad_eulerToMatrix(double alpha, double beta, double gamma,
                          double mat[3][3]);
int cmopad_basis_switcher(enum cmopad_basis_enum in_system,
                          enum cmopad_basis_enum out_system, double r[3][3]);
int cmopad_Eigenvector2PrincipalAxis(enum cmopad_basis_enum coord, double eig, 
                                     double ev[3], double paxis[3]);
int cmopad_basis_transformVector(double v[3],
                                 enum cmopad_basis_enum in_sys,
                                 enum cmopad_basis_enum out_sys);
int cmopad_basis_transformMatrixM6(double *m,
                                   enum cmopad_basis_enum in_sys,
                                   enum cmopad_basis_enum out_sys);
int cmopad_basis_transformMatrixM33(double m[3][3],
                                    enum cmopad_basis_enum in_sys,
                                    enum cmopad_basis_enum out_sys);
void cmopad_strikeDipRake2MT6(double strike, double dip, double rake,
                              double moments[6]);
int cmopad_SetupMT(int nmech, double *mech, 
                   enum cmopad_basis_enum input_basisIn,
                   double Mech_out[3][3]);
int cmopad_findFaultPlanes(int iverb, struct cmopad_struct *src);
int cmopad_MT2PrincipalAxisSystem(int iverb, struct cmopad_struct *src);
int cmopad_standardDecomposition(double Min[3][3], struct cmopad_struct *src);
int cmopad_sortEigenvAscending(bool isabs, double e[3], double ev[3][3]);
void cmopad_uniqueEuler(double *alpha, double *beta, double *gamma);
void cmopad_MatrixToEuler(double rotmat[3][3],
                          double *alpha, double *beta, double *gamma);
void cmopad_printMatrix3x3(bool lfact, double m[3][3]);
double cmopad_momentMagnitudeM6(const double M6[6], int *ierr);
double cmopad_momentMagnitudeM3x3(const double M[3][3], int *ierr);

#ifdef __cplusplus
}
#endif
#endif /* _cmopad_h__ */
