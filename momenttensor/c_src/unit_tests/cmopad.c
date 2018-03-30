#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <lapacke.h>
#include <cblas.h>
#endif
#include "cmopad.h"

/*
Copyright (C) 2010
Lars Krieger & Sebastian Heimann

Contact
lars.krieger@zmaw.de  &  sebastian.heimann@zmaw.de

Modifications - 2015
Converted original Python mopad.py to C - Ben Baker (ISTI)

Contact
benbaker@isti.com

#######################################################################

License:

GNU Lesser General Public License, Version 3

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
*/

#define rad2deg 180.0/M_PI
#define epsilon 1.e-13
static double numpy_mod(double a, double b);
static int __cmopad_argsort3(double *a, int *perm);
static int __cmopad_inv3(double Amat[3][3]);
static double __cmopad_determinant3x3(double A[3][3]);
static double __cmopad_trace3(double M[3][3]);
static int __cmopad_Eigs3x3(double a[3][3], int job, double eigs[3]);
static void __cmopad_swap8(double *a, double *b);
static void array_cross3(const double *__restrict__ u,
                         const double *__restrict__ v,
                         double *__restrict__ n);

//============================================================================//
/*!
 * @brief Generates basis for coordinate transform switch 
 *
 * @param[in] in_system    input system: NED, USE, XYZ, NWU
 * @param[in] out_system   output system: NED, USE, XYZ, NWU
 *
 * @param[out] r           rotation matrix [3 x 3] 
 *
 * @result 0 in indicates success
 *
 */
int cmopad_basis_switcher(enum cmopad_basis_enum in_system, 
                          enum cmopad_basis_enum out_system, double r[3][3])
{
    const char *fcnm = "cmopad_basis_switcher\0";
    double From[9], To[9], r9[9]; 
    double alpha = 1.0;
    double beta  = 0.0;
    int n3 = 3;
    int i, j;
    //------------------------------------------------------------------------//
    //
    // null out to and from
    for (i=0; i<9; i++)
    {
        To[i]   = 0.0;
        From[i] = 0.0;
    }
    // Classify in
    if (in_system == NED)
    {
        // Inverse of identity: From[0] = 1.0; From[4] = 1.0; From[8] = 1.0 is:
        From[0] = 1.0; From[4] = 1.0; From[8] = 1.0;
    }
    else if (in_system == USE)
    {
        From[3] =-1.0; From[7] = 1.0; From[2] =-1.0;  
    }
    else if (in_system == XYZ)
    {
        From[3] = 1.0; From[1] = 1.0; From[8] =-1.0; 
    }
    else if (in_system == NWU)
    {
        From[0] = 1.0; From[4] =-1.0; From[8] =-1.0;
    }
    else
    {
        printf("%s: Invalid in coordinate system %d\n", fcnm, in_system);
        return -1;
    }
    // Classify out
    if (out_system == NED)
    {
        To[0] = 1.0; To[4] = 1.0; To[8] = 1.0;
    }
    else if (out_system == USE)
    {
        // Inverse of: To[3] =-1.0; To[7] = 1.0; To[2] =-1.0; is
        To[1] =-1.0; To[5] = 1.0; To[6] =-1.0;
    }
    else if (out_system == XYZ)
    {
        // Inverse of To[3] = 1.0; To[1] = 1.0; To[8] =-1.0; is:
        To[1] = 1.0; To[3] = 1.0; To[8] =-1.0;
    }
    else if (out_system == NWU)
    {
        //Inverse of: To[0] = 1.0; To[4] =-1.0; To[8] =-1.0; is
        To[0] = 1.0; To[4] =-1.0; To[8] =-1.0;
    }
    else
    {
        printf("%s: Invalid out coordinate system %d\n", fcnm, out_system);
        return -1;
    }  
    // Multiply inv(in)*out
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                n3, n3, n3, alpha, From, n3, To, n3, beta,
                r9, n3);
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            r[i][j] = r9[3*j + i];
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Transform a matrix from in_sys to out_sys.  This is an interface
 *        routine to cmopad_basis__transformMatrixM33
 *        because many matrices are stored as a length 6 array.
 *
 * @param[in,out] m     on input length 6 matrix in in_sys.  the ordering is 
 *                      \f$ {1,2,3,4,5,6} = {11,22,33,12,13,23} \f$.
 *                      on output length 6 matrix in out_sys.  the ordering is
 *                      \f$ {1,2,3,4,5,6} = {11,22,33,12,13,23} \f$.
 *
 * @param[in] in_sys    input system: NED, USE, XYZ, NWU
 * @param[in] out_sys   output system: NED, USE, XYZ, NWU
 *
 * @result 0 indicates success
 *
 */
int cmopad_basis_transformMatrixM6(double *m,
                                   enum cmopad_basis_enum in_sys,
                                   enum cmopad_basis_enum out_sys)
{
    const char *fcnm = "cmopad_basis_transformMatrixM6\0";
    double M33[3][3];
    int ierr;
    //------------------------------------------------------------------------//
    //
    // Convert array matrix to matrix matrix
    M33[0][0] = m[0];             //mxx; mrr
    M33[1][1] = m[1];             //myy; mtt
    M33[2][2] = m[2];             //mzz; mpp
    M33[0][1] = M33[1][0] = m[3]; //mxy; mrt
    M33[0][2] = M33[2][0] = m[4]; //mxz; mrp
    M33[1][2] = M33[2][1] = m[5]; //myz; mtp
    // Switch basis 
    ierr = cmopad_basis_transformMatrixM33(M33, in_sys, out_sys);
    if (ierr != 0)
    {
        printf("%s: Error changing basis\n", fcnm);
        return -1;
    }
    // Convert matrix matrix back to array matrix
    m[0] = M33[0][0]; //mxx; mrr
    m[1] = M33[1][1]; //myy; mtt
    m[2] = M33[2][2]; //mzz; mpp
    m[3] = M33[0][1]; //mxy; mrt
    m[4] = M33[0][2]; //mxz; mrp
    m[5] = M33[1][2]; //myz; mtp
    return 0;
}
//============================================================================//
/*!
 * @brief Transform a symmetric 3x3 matrix from in_sys to out_sys. 
 *
 * @param[in,out] m     on input [3x3] symmetric matrix in in_sys.
 *                      on output [3x3] symmetric matrix in out_sys
 *
 * @param[in] in_sys    input system: NED, USE, XYZ, NWU
 * @param[in] out_sys   output system: NED, USE, XYZ, NWU
 *
 * @result 0 indicates success
 *
 */
int cmopad_basis_transformMatrixM33(double m[3][3],
                                    enum cmopad_basis_enum in_sys,
                                    enum cmopad_basis_enum out_sys)
{
    const char *fcnm = "cmopad_basis_transformMatrixM33\0";
    double r[3][3], t9[9], invR[3][3], invR9[9], r9[9], m9[9];
    double alpha = 1.0;
    double beta = 0.0;
    int n3 = 3;
    int i, j, ierr;
    //------------------------------------------------------------------------//
    //
    // Compute change of basis matrix 
    ierr = cmopad_basis_switcher(in_sys, out_sys, r);
    if (ierr != 0)
    {
        printf("%s: Error generating basis switching matrix\n",fcnm);
        return -1;
    }
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            r9[3*j + i] = r[i][j];
            m9[3*j + i] = m[i][j];
            //invR9[3*j + i] = r[i][j];
            invR[i][j] = r[i][j];
        }
    }
    // compute inv(R)
    ierr = __cmopad_inv3(invR); //double Amat[3][3])
    if (ierr != 0)
    {
        printf("%s: Unlikely error!\n", fcnm);
        return -1;
    }
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            invR9[3*j + i] = invR[i][j];
        }
    }
    // Compute R M inv(R)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                n3, n3, n3, alpha, m9, n3, invR9, n3, beta,
                t9, n3);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                n3, n3, n3, alpha, r9, n3, t9, n3, beta,
                m9, n3);
    // Copy back
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            m[i][j] = m9[3*j + i];
        }
    }
    return 0;
}
//============================================================================//
/*!  
 * @brief Transforms a vector v from basis in_sys to basis out_sys
 *
 * @param[in,out] v     on input length 3 vector in in_sys.  
 *                      on output length 3 vector in out_sys.
 *
 * @param[in] in_sys    input system: NED, USE, XYZ, NWU
 * @param[in] out_sys   output system: NED, USE, XYZ, NWU
 *
 * @result 0 indicates success
 *
 */
int cmopad_basis_transformVector(double v[3],
                                 enum cmopad_basis_enum in_sys,
                                 enum cmopad_basis_enum out_sys)
{
    const char *fcnm = "cmopad_basis_transformVector\0";
    double r[3][3], r9[9], x[3]; 
    double alpha = 1.0;
    double beta = 0.0;
    int n3 = 3;
    int incx = 1;
    int incy = 1;
    int i,j,ierr;
    //------------------------------------------------------------------------//
    //
    // Compute change of basis matrix
    ierr = cmopad_basis_switcher(in_sys, out_sys, r);
    if (ierr != 0)
    {
        printf("%s: Error making basis switcher matrix\n", fcnm);
        return -1;
    }
    for (i=0; i<3; i++)
    {
        x[i] = v[i];
        for (j=0; j<3; j++)
        {
            r9[3*j + i] = r[i][j];
        }
    }
    // Compute v <- R*v 
    cblas_dgemv(CblasColMajor, CblasNoTrans, n3, n3, alpha, r9,
                n3, x, incx, beta, v, incy);
    return 0;
}
//============================================================================//
/*!
 * @brief Cross product of 2 length 3 vectors u and v.  Returns normal vector n.
 * 
 * @param[in] u    first vector in cross product u x v [3]
 * @param[in] v    second vector in cross product u x v [3]
 * 
 * @param[out] n   normal vector n = u x v [3]
 *
 */
static void array_cross3(const double *__restrict__ u,
                         const double *__restrict__ v,
                         double *__restrict__ n)
{
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2]; 
    n[2] = u[0]*v[1] - u[1]*v[0];
    return;
}
//============================================================================//
/*! 
 * @brief Converts the unit eigenvector of length eig in coordinate system
 *        coord to length, azimuth, plunge.  
 *        http://docs.obspy.org/_modules/obspy/imaging/beachball.html#PrincipalAxis
 *        This is an extension of MoPaD
 *
 * @param[in] coord    coordinate system (e.g., NED,USE,XYZ)
 * @param[in] eig      length of eigenvector
 * @param[in] ev       eigenvector in system coordinate system [3]
 * 
 * @param[out] paxis   principal axis store azimuth, plunge, length [3]
 * 
 * @result 0 indicates success
 *
 */
int cmopad_Eigenvector2PrincipalAxis(enum cmopad_basis_enum coord, double eig, 
                                     double ev[3], double paxis[3])
{  
    const char *fcnm = "cmopad_Eigenvector2PrincipalAxis\0";
    double v[3], az, plunge;
    double twopi = 2.0*M_PI;
    double pi180i = 180.0/M_PI;
    int ierr;
    //------------------------------------------------------------------------//
    //
    // Convert eigenvector ev to USE convention 
    cblas_dcopy(3, ev, 1, v, 1);
    ierr = cmopad_basis_transformVector(v, coord, USE);
    if (ierr != 0)
    {
        printf("%s: Error switching basis\n", fcnm);
        return -1; 
    }   
    // Compute azimuth and plunge angles  
    plunge = asin(-v[0]);       //arcsin(-z/r)
    az     = atan2(v[2],-v[1]); //atan(x/y)
    if (plunge <= 0.0)
    {
        plunge =-plunge;
        az = az + M_PI;
    }   
    if (az < 0.0){az = az + twopi;}   //Shift to [0,360]
    if (az > twopi){az = az - twopi;} //Shift to [0,360]
    // Convert to degrees
    plunge = plunge*pi180i;
    az = az*pi180i;
    // Copy back result
    paxis[0] = az;     //First value is azimuth (degrees)
    paxis[1] = plunge; //Second value is plunge (degrees)
    paxis[2] = eig;    //Final value is eigenvalue (distance)
    return 0;
}
//============================================================================//
/*!  
 * @brief Given coordinate system (x,y,z) and rotated system (xs,ys,zs)
 *        the line of nodes is the intersection between the x-y and the xs-ys
 *        planes.
 *
 * @param[in] alpha     is the angle between the z-axis and the zs-axis
 *                      and represents the dip angle (radians)
 * @param[in] beta      is the angle between the x-axis and the line of nodes
 *                      and represents the strike (radians)
 * @param[in] gamma     is the angle between the line of nodes and the xs-axis
 *                      and represents the rake angle (radians)
 *
 * @param[out] mat      rotation matrix built from euler angles [3 x 3]
 *
 * @note Original Python usage for moment tensors:
 *
 * m_unrot = numpy.matrix([[0,0,-1],[0,0,0],[-1,0,0]])
 * rotmat = cmopadEulerToMatrix(dip,strike,-rake)
 * m = rotmat.T * m_unrot * rotmat'''
 *
 */
void cmopad_eulerToMatrix(double alpha, double beta, double gamma,
                          double mat[3][3])
{
    double ca, cb, cg, sa, sb, sg;
    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);

    mat[0][0] = cb*cg-ca*sb*sg; mat[0][1] = sb*cg+ca*cb*sg; mat[0][2] = sa*sg;
    mat[1][0] =-cb*sg-ca*sb*cg; mat[1][1] =-sb*sg+ca*cb*cg; mat[1][2] = sa*cg;
    mat[2][0] = sa*sb;          mat[2][1] =-sa*cb;          mat[2][2] = ca;
    return;
}
//============================================================================//
/*!
 * @brief Computes the moment magnitude from a 6 vector
 *
 * @param[in] M6     moment tensor stored:
 *                    \f$ \{m_{11},m_{22},m_{33},m_{12},m_{13},m_{23} \} \f$
 *                   with the moment tensor terms of units Newton-meters
 *
 * @param[out] ierr  0 indicates success
 *
 * @result moment magnitue
 *
 * @author Ben Baker (ISTI)
 *
 */
double cmopad_momentMagnitudeM6(const double M6[6], int *ierr)
{
    const char *fcnm = "cmopad_momentMagnitudeM6\0";
    double mag;
    const double two_third = 2.0/3.0;
    *ierr = 0;
    mag = M6[0]*M6[0] + M6[1]*M6[1] + M6[2]*M6[2]
        + 2.0*(M6[3]*M6[3] + M6[4]*M6[4] + M6[5]*M6[5]);
    if (mag == 0.0)
    {
        *ierr = 1;
        printf("%s: Error magnitude is zero\n", fcnm);
    }
    else
    {
        mag = two_third*(log10(mag) - 9.1);
    }
    return mag;
}
//============================================================================//
/*!
 * @brief Computes the moment magnitude from a the symmetric moment tensor
 *
 * @param[in] M      symmetric moment tensor with terms given in Newton-meters
 *
 * @param[out] ierr  0 indicates success
 *
 * @result moment magnitue
 *
 * @author Ben Baker (ISTI)
 *
 */
double cmopad_momentMagnitudeM3x3(const double M[3][3], int *ierr)
{
    const char *fcnm = "cmopad_momentMagnitudeM33\0";
    double M6[6], mag;
    M6[0] = M[0][0];
    M6[1] = M[1][1];
    M6[2] = M[2][2];
    M6[3] = M[0][1];
    M6[4] = M[0][2];
    M6[5] = M[1][2];
    mag = cmopad_momentMagnitudeM6(M6, ierr);
    if (*ierr != 0) 
    {
        printf("%s: Error computing moment magnitude\n", fcnm);
    }
    return mag;
}
//============================================================================//
/*!
 * @brief Sets the two angle-triples, describing the faultplanes of the
 *        Double Couple, defined by the eigenvectors P and T of the
 *        moment tensor object.
 *
 *        Define a reference Double Couple with strike = dip =
 *        slip-rake = 0, the moment tensor object's DC is transformed
 *        (rotated) w.r.t. this orientation. The respective rotation
 *        matrix yields the first fault plane angles as the Euler
 *        angles. After flipping the first reference plane by
 *        multiplying the appropriate flip-matrix, one gets the second fault
 *        plane's geometry.
 *
 *        All output angles are in degree
 *
 * @param[in] iverb     controls verbosity for debugging. 0 should be quiet.
 *
 * @param[in,out] src   on input holds plot_clr_order 
 *                      on output holds the strike, dips, and rakes
 *                      describing fault planes fp1 and fp2
 *
 * @result 0 indicates success
 *
 */
int cmopad_findFaultPlanes(int iverb, struct cmopad_struct *src)
{
    const char *fcnm = "cmopad_findFaultPlanes\0";
    int n3 = 3;
    double rot_matrix_fp1[3][3], rot_matrix_fp2[3][3],
           refDC_evecs[9], flip_dc[9], work1[9], work2[9],
           pnt_sorted_EV_matrix[3][3],
           det;
    double sqrt22 = 0.7071067811865476; //sqrt(2)/2 
    double alpha = 1.0, beta = 0.0;
    int i, j;
    //------------------------------------------------------------------------//
    //
    // reference Double Couple (in NED) basis w/ (strike,dip,slip) = (0,0,0)
    refDC_evecs[0] =-sqrt22; refDC_evecs[3] = 0.0; refDC_evecs[6] =-sqrt22; 
    refDC_evecs[1] = 0.0;    refDC_evecs[4] = 1.0; refDC_evecs[7] = 0.0;
    refDC_evecs[2] =-sqrt22; refDC_evecs[5] = 0.0; refDC_evecs[8] = sqrt22;

    //refDC_evals[0] =-1.0; 
    //refDC_evals[1] = 0.0;
    //refDC_evals[2] = 1.0;
    // Matrix which is turning from one fault plane to the other
    flip_dc[0] = 0.0; flip_dc[3] = 0.0; flip_dc[6] =-1.0; 
    flip_dc[1] = 0.0; flip_dc[4] =-1.0; flip_dc[7] = 0.0;
    flip_dc[2] =-1.0; flip_dc[5] = 0.0; flip_dc[8] = 0.0;
    // Euler-tools need matrices of EV sorted in PNT (pressure,null,tension):
    for (i=0; i<3; i++)
    {
       for (j=0; j<3; j++)
       {
          pnt_sorted_EV_matrix[i][j] = src->rotation_matrix[i][j];
       }
    } 
    // Re-sort only necessary if abs(p) <= abs(t)
    if (src->plot_clr_order < 0)
    {
        for (i=0; i<3; i++)
        {
            pnt_sorted_EV_matrix[i][0] = src->rotation_matrix[i][2];
            pnt_sorted_EV_matrix[i][2] = src->rotation_matrix[i][0];
        }
    }
    // Rotation matrix, describing the rotation of the eigenvector
    // system of the input moment tensor into the eigenvector
    // system of the reference Double Couple
    if (iverb > 1)
    {
        printf("%s: Pressure, null, tension eigenvectors:\n",fcnm);
    }
    for (j=0; j<3; j++)
    {
       for (i=0; i<3; i++)
       {
          work1[j*n3 + i] = pnt_sorted_EV_matrix[i][j];
       }
       if (iverb > 1)
       {
           printf(
           "%8.5f %8.5f %8.5f\n",pnt_sorted_EV_matrix[j][0],
                                 pnt_sorted_EV_matrix[j][1],
                                 pnt_sorted_EV_matrix[j][2]);
       }
    }
    if (iverb > 1)
    {
        printf("%s: Reference double couple eigenvectors:\n",fcnm);
        for (i=0; i<3; i++)
        {
            printf(
            "%8.5f %8.5f %8.5f\n",refDC_evecs[i],
                                  refDC_evecs[3+i],refDC_evecs[6+i]);
        }
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n3, n3, n3, 
                alpha, work1, n3, refDC_evecs, n3, beta,
                work2, n3);
    for (j=0; j<3; j++)
    {
        for (i=0; i<3; i++)
        {
            rot_matrix_fp1[j][i] = work2[j*n3 + i]; //Copy transpose
        }
    }
    // check if rotation has correct orientation
    det = __cmopad_determinant3x3(rot_matrix_fp1);
    if (det < 0.0)
    {
       for (i=0; i<3; i++)
          {
          for (j=0; j<3; j++)
          {
              rot_matrix_fp1[i][j] =-1.0*rot_matrix_fp1[i][j];
          }
       }
    } 
    for (j=0; j<3; j++)
    {
        for (i=0; i<3; i++)
        {
            work1[j*3 + i] = rot_matrix_fp1[i][j];
        }
    }
    // adding a rotation into the (ambiguous) system of the 2nd fault plane
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3, n3,
                alpha, flip_dc, n3, work1, n3, beta, 
                work2, n3);
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            rot_matrix_fp2[i][j] = work2[j*3 + i];
        }
    }
    // Calculate strike, dip, rake of fault planes
    cmopad_findStrikeDipRake(rot_matrix_fp1, 
                             &src->fp1[1], &src->fp1[0], &src->fp1[2]);
    cmopad_findStrikeDipRake(rot_matrix_fp2,
                             &src->fp2[1], &src->fp2[0], &src->fp2[2]);
    // For convenience if choose the first fault plane to be the one with
    // the smaller strike
    if (src->fp2[0] < src->fp1[0])
    {
        for (i=0; i<3; i++)
        {
            __cmopad_swap8(&src->fp1[i], &src->fp2[i]);
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Returns dip, strike, and slip(rake) angles in degrees describing the 
 *        fault plane.
 *
 * @param[in] rot_mat     matrix of eignevectors for fault plane [3 x 3]
 *
 * @param[out] alpha      dip angle (degrees)
 * @param[out] beta       strike angle (degrees)
 * @param[out] gamma      rake angle (degrees)
 *
 */
void cmopad_findStrikeDipRake(double rot_mat[3][3],
                              double *alpha, double *beta, double *gamma)
{
    double a = 0.0, b = 0.0, g = 0.0;
    cmopad_MatrixToEuler(rot_mat, &a, &b, &g);
    *alpha = a*rad2deg;
    *beta  = b*rad2deg;
    *gamma =-g*rad2deg;
    return;
}
//============================================================================//
/*!
 * @brief Returns three Euler angles alpha, beta, and gamma (radians) from a
 *        rotation matrix.
 *
 * @param[in] rotmat    rotation matrix of which to decompose into euler angles
 *                      [3 x 3]
 *
 * @param[out] alpha    euler angle representing dip (radians)
 * @param[out] beta     euler angle representing strike (radians)
 * @param[out] gamma    euler angle representing -rake (radians)
 *
 */
void cmopad_MatrixToEuler(double rotmat[3][3],
                          double *alpha, double *beta, double *gamma)
{
    const int n3 = 3;
    int inc = 1;
    double ex[3], ez[3], exs[3], ezs[3], enodes[3], enodess[3], 
           xnorm, cos_alpha, a,b,g,a1,a2;
    double twopi = 2.0*M_PI;
    int i,j;
    //------------------------------------------------------------------------//
    //
    // transpose matrix vector multiply
    ex[0] = 1.0; ex[1] = 0.0; ex[2] = 0.0;
    ez[0] = 0.0; ez[1] = 0.0; ez[2] = 1.0;
    for (i=0; i<3; i++)
    {
        exs[i] = 0.0; ezs[i] = 0.0;
        for (j=0; j<3; j++)
        {
            exs[i] = exs[i] + rotmat[j][i]*ex[j];
            ezs[i] = ezs[i] + rotmat[j][i]*ez[j];
        }
    }
    // Cross product
    array_cross3(ez, ezs, enodes);
    xnorm = sqrt(pow(enodes[0],2) + pow(enodes[1],2) + pow(enodes[2],2));
    if (xnorm < 1.e-10)
    {
        for (i=0; i<3; i++)
        {
            enodes[i] = exs[i];
        }
    }
    for (i=0; i<3; i++)
    {
        enodess[i] = 0.0;
        for (j=0; j<3; j++)
        {
            enodess[i] = enodess[i] + rotmat[i][j]*enodes[j];
        }
    }
    // Multiply
    cos_alpha = cblas_ddot(n3, ez, inc, ezs, inc); 
    cos_alpha = fmin( 1.0, cos_alpha);
    cos_alpha = fmax(-1.0, cos_alpha);
    a = acos(cos_alpha); 
    a1 = atan2( enodes[1],  enodes[0]);
    a2 =-atan2(enodess[1], enodess[0]);
    b = numpy_mod(a1, twopi);
    g = numpy_mod(a2, twopi);
    // Compute unique Euler angles
    cmopad_uniqueEuler(&a, &b, &g);
    *alpha = a; 
    *beta = b; 
    *gamma = g;
    return;
}
//============================================================================//
/*!
 * @brief Clone of the numpy mod(a,b) function
 *
 * @result numpy.mod(a, b)
 *
 */
static double numpy_mod(double a, double b)
{
    if (b == 0.0){return a;}
    if (a >= 0.0)
    {
        if (b > 0.0)
        {
            return a - floor(a/b)*b;
        }
        else
        {
            return a - fabs(floor(a/b))*fabs(b);
        } 
    }
    else
    {
        if (b > 0.0)
        {
            return a + fabs(floor(a/b))*fabs(b);
        }
        else
        {
            return a + floor(a/b)*fabs(b);
        }
    }
}
//============================================================================//
/*!  
 * Read in Matrix M and set up eigenvalues (EW) and eigenvectors
 * (EV) for setting up the principal axis system.
 *
 * The internal convention is the 'HNS'-system: H is the
 * eigenvector for the smallest absolute eigenvalue, S is the
 * eigenvector for the largest absolute eigenvalue, N is the null
 * axis.
 *
 * Naming due to the geometry: a CLVD is
 * Symmetric to the S-axis,
 * Null-axis is common sense, and the third (auxiliary) axis
 * Helps to construct the RA3.
 *
 * Additionally builds matrix for basis transformation back to NED system.
 *
 * The eigensystem setup defines the colouring order for a later
 * plotting in the BeachBall class. This order is set by the
 * `_plot_clr_order' attribute.
 *
 * @param[in] iverb     controls verbosity.  0 is quiet.
 *
 * @param[in,out] src   on input holds the input moment tensor 
 *                      on output has the eigenvalues and fault planes
 *
 */
int cmopad_MT2PrincipalAxisSystem(int iverb, struct cmopad_struct *src)
{
    const char *fcnm = "cmopad_MT2PrincipalAxisSystem\0";
    int ival = 0;
    double M[3][3], M_devi[3][3], EV_devi[3][3], EV[3][3], EW[3], 
           EW_devi[3], EV1[3], EV2[3], EV3[3], EVh[3], EVs[3], EVn[3],  
           EWs, EWh, EWn,
           EW1, EW2, EW3, EW1_devi, EW2_devi, EW3_devi, 
           trace_M, xsum;
    int EW_order[3], symmetry_around_tension, clr, i, j, ierr;
    //------------------------------------------------------------------------//
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            M[i][j]      = src->M[i][j];
            M_devi[i][j] = src->M_devi[i][j];
            EV[i][j] = M[i][j];
            EV_devi[i][j] = M_devi[i][j];
        }
    }
    ierr = __cmopad_Eigs3x3(EV_devi, ival, EW_devi);
    if (ierr != 0)
    {
        printf("%s: Error computing eigenvalues (a)!", fcnm);
        return -1; 
    } 
    ierr = __cmopad_argsort3(EW_devi, EW_order);
    // Removed if statement, always true in python code
    trace_M = __cmopad_trace3(M); // M[0][0] + M[1][1] + M[2][2];
    if (trace_M < epsilon){trace_M = 0.0;}
    ierr = __cmopad_Eigs3x3(EV, ival, EW);
    if (ierr != 0)
    {
        printf("%s: Error computing eigenvalues (b)!", fcnm);
        return -1;
    }
    for (i=0; i<3; i++)
    {
        if (fabs(EW[i]) < epsilon){EW[i] = 0.0;}
    }
    //trace_M_devi = cmopadTrace(3,M_devi);

    // Save eigenvalues/eigenvectors for deviatoric and full moment tensor
    EW1_devi = EW_devi[EW_order[0]];
    EW2_devi = EW_devi[EW_order[1]];
    EW3_devi = EW_devi[EW_order[2]];
    if (fabs(EW2_devi) < epsilon){EW2_devi = 0.0;}

    EW1 = EW[EW_order[0]];
    EW2 = EW[EW_order[1]];
    EW3 = EW[EW_order[2]];
    for (i=0; i<3; i++)
    {
        //EV1_devi[i] = EV_devi[i][EW_order[0]];
        //EV2_devi[i] = EV_devi[i][EW_order[1]]; 
        //EV3_devi[i] = EV_devi[i][EW_order[2]]; 
        EV1[i] = EV[i][EW_order[0]];
        EV2[i] = EV[i][EW_order[1]];
        EV3[i] = EV[i][EW_order[2]]; 
    }

    // Classify tension axis symmetry
    symmetry_around_tension = 1;
    clr = 1;
    xsum = EW1 + EW2 + EW3;
    // implosion
    if (EW1 < 0.0 && EW2 < 0.0 && EW3 < 0.0)
    {
        if (iverb > 2){printf("%s: Implosion\n", fcnm);}
        symmetry_around_tension = 0;
        clr = 1;
    }
    // explosion
    else if (EW1 > 0.0 && EW2 > 0.0 && EW3 > 0.0)
    {
        if (iverb > 2){printf("%s: Explosion\n", fcnm);}
        symmetry_around_tension = 1;
        if (fabs(EW1_devi) > fabs(EW3_devi))
        {
            symmetry_around_tension = 0;
        }
        clr = -1;
    }
    // net-implosion
    else if (EW2 < 0.0 && xsum < 0.0)
    {
        if (iverb > 2){printf("%s: Net-implosion (1)\n", fcnm);}
        if (fabs(EW1_devi) < fabs(EW3_devi))
        {
            symmetry_around_tension = 1;
            clr = 1; 
        }
        else
        {
            symmetry_around_tension = 1;
            clr = 1;
        }
    }
    // net-implosion
    else if (EW2_devi >= 0.0 && xsum < 0.0)
    {
        if (iverb > 2){printf("%s: Net-impolosion (2)\n", fcnm);}
        symmetry_around_tension = 0;
        clr = -1;
        if (fabs(EW1_devi) < fabs(EW3_devi))
        {
            symmetry_around_tension = 1;
            clr = 1;
        }
    }
    // net-explosion
    else if (EW2_devi < 0.0 && xsum > 0.0)
    {
        if (iverb > 2){printf("%s: Net-explosion (1)\n", fcnm);}
        symmetry_around_tension = 1;
        clr = 1;
        if (fabs(EW1_devi) > fabs(EW3_devi))
        {
            symmetry_around_tension = 0;
            clr = -1;
        }
    }
    // net-explosion
    else if ( EW2_devi >= 0.0 && xsum > 0.0)
    {
        if (iverb > 2){printf("%s: Net-explosion (2)\n", fcnm);}
        symmetry_around_tension = 0;
        clr = -1;
    }
    else
    {
        // If the trace is 0 to machine epsilon then purely deviatoric
        if (trace_M/fmax(fabs(EW1_devi),fabs(EW3_devi)) > 1.e-10){ //!= 0.0){
            printf("%s: Warning should not be here %f %f %f\n",
                      fcnm, EW2_devi, trace_M, xsum);
        }
    }
    if (iverb > 1)
    {
        printf("%s: sym/clr %d %d\n",
                  fcnm, symmetry_around_tension, clr);
    }
    // Deviatoric classification
    if (fabs(EW1_devi) < fabs(EW3_devi))
    {
        symmetry_around_tension = 1;
        clr = 1;
    }
    if (fabs(EW1_devi) >= fabs(EW3_devi))
    {
        symmetry_around_tension = 0;
        clr = -1;
    }
    //explosive src with all neg evals? - TODO: email lars.  an edge
    //case is (-1,-1,-1, 0,0,0) which has trace_M == 0.0 and is a pure
    //implosion.  in this case there is a breakdown of the code because
    //trace_M == 0.0
    if (EW3 < 0.0 && trace_M > 0.0)
    {
        printf("%s: Critical error!\n", fcnm);
        return -1;
    }

    // Classify pure deviatoric - pure deviatoric
    if (trace_M == 0.0)
    {
        if (iverb > 1){printf("%s: Pure-deviatoric\n", fcnm);}
        // pure shear
        if (EW2 == 0.0)
        {
            if (iverb > 1){printf("%s: Pure-shear\n", fcnm);}
            symmetry_around_tension = 1;
            clr = 1;
        }
        else if (2.0*fabs(EW2) == fabs(EW1) || 2.0*fabs(EW2) == fabs(EW3))
        {
            if (iverb > 1){printf("%s: Pure CLVD\n", fcnm);}
            // CLVD symmetry around tension
            if (fabs(EW1) < EW3)
            {
                if (iverb > 2)
                {
                    printf("%s: CLVD symmetry around tension\n", fcnm);
                }
                symmetry_around_tension = 1;
                clr = 1;
            }
            // CLVD symmetry around pressure
            else
            {
                if (iverb > 2)
                {
                    printf("%s: CLVD symmetry around pressure\n", fcnm);
                }
                symmetry_around_tension = 0;
                clr = -1;
            }
        }
        // mix of DC and CLVD
        else
        {
            if (iverb > 1)
            {
                printf("%s: Mix of DC and CLVD\n", fcnm);
            }
            // symmetry around tension
            if (fabs(EW1) < EW3)
            {
                if (iverb > 2)
                {
                    printf("%s: Symmetry around tension\n", fcnm);
                }
                symmetry_around_tension = 1;
                clr = 1;
            }
            // symmetry around pressure
            else
            {
                if (iverb > 2)
                {
                    printf("%s: Symmetry around pressure\n", fcnm);
                }
                symmetry_around_tension = 0;
                clr = -1;
            }
        }
    }
    if (iverb > 0)
    {
        printf("%s: Final symmetry %d %d\n",
                  fcnm, symmetry_around_tension, clr);
    }
    // Define order of eigenvectors and values according to symmetry axis
    if (symmetry_around_tension == 1)
    {
        EWs = EW3;
        EWh = EW1; 
        EWn = EW2;
        for (i=0; i<3; i++)
        {
            EVs[i] = EV3[i];
            EVh[i] = EV1[i];
            EVn[i] = EV2[i];
        }
    }
    else
    {
        EWs = EW1;
        EWh = EW3;
        EWn = EW2;
        for (i=0; i<3; i++)
        {
            EVs[i] = EV1[i];
            EVh[i] = EV3[i];
            EVn[i] = EV2[i]; 
        }
    }
    // Order of eigenvector's basis: (H,N,S)
    for (i=0; i<3; i++)
    {
        // Matrix for basis transform
        src->rotation_matrix[i][0] = EVh[i]; 
        src->rotation_matrix[i][1] = EVn[i];
        src->rotation_matrix[i][2] = EVs[i]; 
        // Save eigenvectors 
        src->eigenvectors[i][0] = EVh[i];
        src->eigenvectors[i][1] = EVn[i];
        src->eigenvectors[i][2] = EVs[i];
        // Principal axis
        src->p_axis[i] = EV1[i];    //Smallest (most negative) eval/evec pair
        src->null_axis[i] = EVn[i]; //Intermediate eval/evec pair
        src->t_axis[i] = EV3[i];    //Largest (most positive) eval/evec pair 
        /* TODO the above is definitionally correct but the following works. 
           What's going on here?  Did someone switch conventions on me? */
        //src->t_axis[i] = EV1[i];
        //src->null_axis[i] = EVn[i];
        //src->p_axis[i] = EV3[i];
    }
    src->eig_pnt[0] = EW1_devi; // Smallest - most negative eigenvalue
    src->eig_pnt[1] = EW2_devi;
    src->eig_pnt[2] = EW3_devi; // Largest eigenvalue
    // TODO - these are reversed to match above t_axis, null_axis, and p_axis
    //printf("%f %f\n", EW1_devi, EW3_devi);
    //src->eig_pnt[0] = EW3_devi; // Make pressure most positive eigenvalue?
    //src->eig_pnt[1] = EW2_devi;
    //src->eig_pnt[2] = EW1_devi; // Make tension most negative eigenvalue
    // Eigenvalues corresponding to eigenvectors (H,N,S) basis
    src->eigenvalues[0] = EWh;
    src->eigenvalues[1] = EWn;
    src->eigenvalues[2] = EWs;

    // Important for plot in BeachBall class (also used in FindFaultPlanes)
    src->plot_clr_order = clr;

    ierr = cmopad_findFaultPlanes(iverb, src);
    if (ierr != 0)
    {
        printf("%s: Error computing fault planes!\n", fcnm);
    }
    return ierr;
}
//============================================================================//
/*!  
 * @brief Prints a 3x3 matrix in a Mopad-esque style.
 *
 * @param[in] lfact   if true then use a normalization factor
 * @param[in] m       the 3x3 matrix to write to standard out
 * 
 */
void cmopad_printMatrix3x3(bool lfact, double m[3][3])
{
    const char *fcnm = "cmopad_printMatrix3x3\0";
    double norm_fact, m11,m22,m33,m12,m13,m23;
    int i,j;
    //------------------------------------------------------------------------//
    norm_fact = 1.0;
    if (lfact)
    {
        norm_fact = 0.0;
        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                //norm_fact = fmax(norm_fact,fabs(mathRound(m[i][j],5)));
                norm_fact = fabs(m[i][j]);
            }
        }
        if (norm_fact == 0.0)
        {
            printf("%s: Warning norm_fact = 0.0!\n", fcnm);
            return;
        }
    }
    m11 = m[0][0]/norm_fact;
    m22 = m[1][1]/norm_fact;
    m33 = m[2][2]/norm_fact;
    m12 = m[0][1]/norm_fact; 
    m13 = m[0][2]/norm_fact;
    m23 = m[1][2]/norm_fact; 
    if (lfact)
    {
        printf(" [ %5.2f %5.2f %5.2f ] \n",m11,m12,m13);
        printf(" [ %5.2f %5.2f %5.2f ] x %f\n",
                          m12,m22,m23,norm_fact);
        printf(" [ %5.2f %5.2f %5.2f ] \n",m13,m23,m33);
    }
    else
    {
        printf(" [ %5.2f %5.2f %5.2f ]\n",m11,m12,m13);
        printf(" [ %5.2f %5.2f %5.2f ]\n",m12,m22,m23);
        printf(" [ %5.2f %5.2f %5.2f ]\n",m13,m23,m33);
    }
    return;
}
//============================================================================//
/*!
 * @brief Brings the provided mechanism into symmetric 3x3 matrix form.
 *
 *        The source mechanism may be provided in different forms:
 *
 *        -- as 3x3 matrix - symmetry is checked - one basis 
 *           system has to be chosen, or NED as default is taken
 *        -- as 3-element tuple or array - interpreted as 
 *           strike, dip, slip-rake angles in degree
 *        -- as 4-element tuple or array - interpreted as strike, dip, 
 *           slip-rake angles in degree + seismic scalar moment in Nm
 *        -- as 6-element tuple or array - interpreted as the 6 independent 
 *           entries of the moment tensor
 *        -- as 7-element tuple or array - interpreted as the 6 independent 
 *           entries of the moment tensor + seismic scalar moment in Nm
 *        -- as 9-element tuple or array - interpreted as the 9 entries 
 *           of the moment tensor - checked for symmetry
 *        -- as a nesting of one of the upper types (e.g. a list of 
 *           n-tuples) - first element of outer nesting is taken
 *
 */
int cmopad_SetupMT(int nmech, double *mech, 
                   enum cmopad_basis_enum input_basisIn,
                   double Mech_out[3][3])
{
    const char *fcnm = "cmopad_SetupMT\0";
    enum cmopad_basis_enum input_basis;
    double rotmat[9], mtemp1[9], mtemp2[9], m_unrot[9], scalar_moment,  
           strike, dip, rake;
    int i, j, ierr;
    double pi180 = M_PI/180.0;
    double alpha = 1.0;
    double beta = 0.0;
    int n3 = 3;
    //------------------------------------------------------------------------//
    //
    // All 9 elements are specified
    input_basis = input_basisIn;
    if (nmech == 9)
    {
        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                Mech_out[i][j] = mech[3*j + i];
            }
        }
    }
    // Mechanism given as 6 or 7 element array
    else if (nmech == 6 || nmech == 7)
    {
        scalar_moment = 1.0;
        if (nmech == 7){scalar_moment = mech[6];}
        Mech_out[0][0] = scalar_moment*mech[0];
        Mech_out[1][1] = scalar_moment*mech[1];
        Mech_out[2][2] = scalar_moment*mech[2];
        Mech_out[0][1] = Mech_out[1][0] = scalar_moment*mech[3];
        Mech_out[0][2] = Mech_out[2][0] = scalar_moment*mech[4];
        Mech_out[1][2] = Mech_out[2][1] = scalar_moment*mech[5];
    }
    // Strike dip rake; conventions from Jost and Herrmann; NED-basis
    else if (nmech == 3 || nmech == 4)
    {
        scalar_moment = 1.0;
        if (nmech == 4){scalar_moment = mech[3];}
        // Set matrix as function of Euler angles
        strike = mech[0]*pi180;
        dip    = mech[1]*pi180;
        rake   = mech[2]*pi180;
        cmopad_eulerToMatrix(dip, strike, -rake, Mech_out);
        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                m_unrot[3*j + i] = 0.0;
                rotmat[3*j + i] = Mech_out[i][j];
                mtemp1[3*j + i] = 0.0;
                mtemp2[3*j + i] = 0.0;
            }
        }
        m_unrot[6] =-1.0;
        m_unrot[2] =-1.0;
        // Multiply trans(Rotmat)*m_unrot*Rotmat 
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n3, n3, n3,
                    alpha, rotmat, n3, m_unrot, n3, beta,
                    mtemp1, n3);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3, n3,
                    alpha, rotmat, n3, mtemp1, n3, beta,
                    mtemp2, n3);
        // Copy back and include scalar moment
        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                Mech_out[i][j] = scalar_moment*mtemp2[3*j + i]; 
            }
        }
        // Assure right basis system
        input_basis = NED;
    }
    // Unknown mechanism
    else
    {
        printf("%s: Invalid mechanism!\n", fcnm);
        return -1;
    }
    // Rotate into NED basis for local processing
    ierr = cmopad_basis_transformMatrixM33(Mech_out, input_basis, NED);
    if (ierr != 0)
    {
        printf("%s: Error transforming matrix\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!  
 * @brief Sorts eigenvalues/eigenvectors of a matrix into ascending order.
 *        Optionally, can sort based on the absolute magnitude of eigenvalues
 *
 * @param[in] isabs     if true then sort on absolute value of eigenvalues
 * 
 * @param[in,out] e     eigenvalues sorted in ascending order [3]
 * @param[in,out] ev    corresponding eigenvectors sorted along with eigenvalues
 *                      [3 x 3]
 *
 * @result 0 indicates success
 *
 */
int cmopad_sortEigenvAscending(bool isabs, double e[3], double ev[3][3])
{
    const char *fcnm = "cmopad_sortEigenvAscending\0";
    const int n = 3;
    double evwork[3][3], ework[3], esave[3];
    int iperm[3], i, ierr, j;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        esave[i] = e[i];
        if (isabs)
        {
            ework[i] = fabs(e[i]);
        }
        else
        {
            ework[i] = e[i]; 
        }
        for (j=0; j<n; j++)
        {
            evwork[i][j] = ev[i][j];
        }
    }
    // Sort
    ierr = __cmopad_argsort3(ework, iperm);
    if (ierr != 0)
    {
        printf("%s: Error in permutation sort\n", fcnm);
        return -1;
    }
    // Copy back into ascending order
    for (i=0; i<n; i++)
    {
        e[i] = esave[iperm[i]];
        for (j=0; j<n; j++)
        {
            ev[i][j] = evwork[i][iperm[j]];
        }
    }
    return ierr;
}
//============================================================================//
/*!  
 * @brief Decomposes the input moment tensor Min in the North, East, Down
 *        coordinates into an isotropic, compensated linear vector dipole, 
 *        and pure double couple.  
 *
 * @param[in] Min   3 x 3 moment tensor in NED coordinates to decompose
 *                  into a pure isotropic, CLVD, and DC  
 *
 * @param[out] src  internal mopad structure.  this routine will compute:
 *                  M:      copy of the original moment tensor
 *                  M_iso:  diagonal isotropic moment tensor
 *                  M_devi: deviatoric moment tensor (M - M_iso) 
 *                  M_CLVD: CLVD contribution to M_devi
 *                  M_DC:   DC contribution to M_devi 
 *                  DC_percentage:   Percent (0,100) double couple
 *                  CLVD_percentage: Percent (0,100) CLVD
 *                  ISO_Percetnage:  Percent (0,100) isotropic
 *                  seismic_moment:  seismic moment (Nm) Bowers & Hudson 
 *                                   1999
 *                  moment_magnitude: moment magnitude Hanks and Kanamori
 *                                    1979 
 * 
 * @result 0 indicates success
 *
 */
int cmopad_standardDecomposition(double Min[3][3], struct cmopad_struct *src)
{
    const char *fcnm = "cmopad_standardDecomposition\0";
    double M[3][3], M_devi[3][3], M_iso[3][3],   
           M_DC[3][3], M_CLVD[3][3], eigenvtot[3][3], eigenvdevi[3][3],
           a3out, a2out, a1out, eigenwtot[3], eigenwdevi[3], 
           a1[3], a2[3], a3[3];
    double M0, M0_iso, M0_devi, F, M_DC_percentage, M_iso_percentage, tracem;
    int lsub, i,j, ierr;
    int ival = 0; //Calculate eigenvalues only (1)   
    // Copy 
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            M[i][j] = Min[i][j];
        }
    }
    // Isotropic part
    tracem = __cmopad_trace3(M); //M[0][0] + M[1][1] + M[2][2];
    M_iso[0][0] = M_iso[0][1] = M_iso[0][2] = 0.0;
    M_iso[1][0] = M_iso[1][1] = M_iso[1][2] = 0.0;
    M_iso[2][0] = M_iso[2][1] = M_iso[2][2] = 0.0;
    M_iso[0][0] = M_iso[1][1] = M_iso[2][2] = (1.0/3.0)*tracem;
    M0_iso = fabs(1.0/3.0*tracem);

    // Deviatoric part
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            M_devi[i][j] = M[i][j] - M_iso[i][j]; 
        }
    }
    // Copy to structure
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            src->M[i][j] = M[i][j];
            src->M_iso[i][j]  = M_iso[i][j];
            src->M_devi[i][j] = M_devi[i][j];
        }
    }
 
    // Eigenvalues and eigenvectors
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            eigenvtot[i][j] = M_devi[i][j]; //I dont understand this
            eigenvdevi[i][j]= M_devi[i][j];
        }
    }
    ierr = __cmopad_Eigs3x3(eigenvtot, ival, eigenwtot);
    if (ierr != 0)
    {
        printf("%s: Error computing eigs (a)!\n",fcnm);
        return -1;
    } 
    ierr = __cmopad_Eigs3x3(eigenvdevi, ival, eigenwdevi);
    if (ierr != 0)
    {
        printf("%s: Error computing eigs (b)!\n",fcnm);
        return -1;
    }

    // Eigenvalues in ascending order
    ierr = ierr + cmopad_sortEigenvAscending(true, eigenwtot , eigenvtot);
    ierr = ierr + cmopad_sortEigenvAscending(true, eigenwdevi, eigenvdevi);
    if (ierr != 0)
    {
        printf("%s: Error sorting eigenvalues/eigenvectors\n", fcnm);
        return -1;
    }

    /* TODO: Email Lars Krieger; confirm this is a mistake
           : Jost and Herrmann Eqn 19 
    M0_devi =-1.0;
    for (i=0; i<3; i++){
        M0_devi = fmax(M0_devi, fabs(eigenwdevi[i]));
    }
    */
    // Compute mangnitude (Jost and Herrmann Eqn 19 noting the eigenvalues are
    // in ascending order)
    M0_devi = 0.5*(fabs(eigenwdevi[1]) + fabs(eigenwdevi[2])); 

    // Named according to Jost and Herrmann
    for (i=0; i<3; i++)
    {
        a1[i] = eigenvtot[i][0];
        a2[i] = eigenvtot[i][1];
        a3[i] = eigenvtot[i][2];
    }

    // If only isotropic part exists
    if (M0_devi < epsilon)
    {
        F = 0.5;
    }
    else
    {
        F =-eigenwdevi[0]/eigenwdevi[2]; 
    }

    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            a3out = a3[i]*a3[j];
            a2out = a2[i]*a2[j];
            a1out = a1[i]*a1[j];
            // Second term of eqn 37 Jost and Herrmann
            M_DC[i][j] = eigenwdevi[2]*(1.0 - 2.0*F)*(a3out - a2out);
            // Remainder of deviatoric is CLVD (could do 3rd term of eqn 37) 
            lsub = 1;
            if (lsub == 1)
            {
                M_CLVD[i][j] = M_devi[i][j] - M_DC[i][j];
            }
            else
            {
                M_CLVD[i][j] = eigenwdevi[2]*F
                              *(2.0*a3out - a2out - a1out);
            }
        }
    }

    //according to Bowers & Hudson:
    M0 = M0_iso + M0_devi;

    // Really shouldn't be rounding 
    /* 
    M_iso_percentage     = mathRound(M0_iso/M0*100.0,6);
    //cmopad.iso_percentage = M_iso_percentage;
    M_DC_percentage = ( mathRound( (1.0 - 2.0*fabs(F))
                                  *(1.0 - M_iso_percentage/100.0)*100,6));
    */
    M_iso_percentage = M0_iso/M0*100.0;
    M_DC_percentage = (1.0 - 2.0*fabs(F))*(1.0 - M_iso_percentage/100.0)*100.0;
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            src->DC[i][j]   =  M_DC[i][j];
            src->CLVD[i][j] =  M_CLVD[i][j];
        }
    }
    src->DC_percentage   = M_DC_percentage;
    src->ISO_percentage  = M_iso_percentage;
    src->DEV_percentage  = 100.0 - src->ISO_percentage;
    src->CLVD_percentage = 100.0 - src->ISO_percentage - src->DC_percentage;
    ////src->seismic_moment   = sqrt(1./2.0*(pow(eigenw[0],2) + pow(eigenw[0],2) + pow(eigenw[0],2))); 
    src->seismic_moment   = M0;
    src->moment_magnitude = log10(src->seismic_moment*1.0e7)/1.5 - 16.1/1.5;

    return 0;
}
//============================================================================//
/*!  
 * @brief Return 6-tuple containing entries of M, calculated from fault plane 
 *        angles (defined as in Jost&Herman), given in degrees.
 *        Basis for output is NED (= X,Y,Z)
 *
 * @param[in] strike    angle clockwise between north and plane (in [0,360])
 * @param[in] dip       angle between surface and dipping plane (in [0,90]) 
 *                      0 = horizontal, 90 = vertical
 * @param[in] rake      angle on the rupture plane between strike vector 
 *                      and actual movement (defined mathematically 
 *                      counter-clockwise rotation is positive)
 *
 * @param[out] moments  M = {M_nn, M_ee, M_dd, M_ne, M_nd, M_ed}
 *
 */
void cmopad_strikeDipRake2MT6(double strike, double dip, double rake,
                              double moments[6])
{
    double M1,M2,M3,M4,M5,M6;
    double S_rad = strike/rad2deg;
    double D_rad = dip/rad2deg;
    double R_rad = rake/rad2deg;
    //-------------------------------------------------------------------------//
    if (fabs(S_rad) < epsilon){S_rad = 0.0;}
    if (fabs(D_rad) < epsilon){D_rad = 0.0;}
    if (fabs(R_rad) < epsilon){R_rad = 0.0;}

    M1 =-( sin(D_rad)*cos(R_rad)*sin(2.0*S_rad) 
         + sin(2.0*D_rad)*sin(R_rad)*pow(sin(S_rad),2) );
    M2 = ( sin(D_rad)*cos(R_rad)*sin(2.0*S_rad) 
         - sin(2.0*D_rad)*sin(R_rad)*pow(cos(S_rad),2) );
    M3 = sin(2.0*D_rad)*sin(R_rad);
    M4 = sin(D_rad)*cos(R_rad)*cos(2.0*S_rad) 
       + 0.5*sin(2.0*D_rad)*sin(R_rad)*sin(2.0*S_rad);
    M5 =-( cos(D_rad)*cos(R_rad)*cos(S_rad) 
         + cos(2.0*D_rad)*sin(R_rad)*sin(S_rad) );
    M6 =-( cos(D_rad)*cos(R_rad)*sin(S_rad)   
         - cos(2.0*D_rad)*sin(R_rad)*cos(S_rad));

    moments[0] = M1; 
    moments[1] = M2; 
    moments[2] = M3; 
    moments[3] = M4; 
    moments[4] = M5; 
    moments[5] = M6;
    return;
}
//============================================================================//
/*!
 * @brief Uniquify euler angle triplet.
 *
 * Puts euler angles into ranges compatible with (dip,strike,-rake) in
 * seismology:
 *
 * @param[in,out] alpha    dip angle to put into range [0, pi/2]
 * @param[in,out] beta     strike angle to put into range [0, 2*pi)
 * @param[in,out] gamma    -rake angle to put into range [-pi, pi)
 *
 * @note If alpha is near to zero, beta is replaced by beta+gamma and gamma is
 *       set to zero, to prevent that additional ambiguity.
 *
 *       If alpha is near to pi/2, beta is put into the range [0,pi).
 *
 */
void cmopad_uniqueEuler(double *alpha, double *beta, double *gamma)
{
    const char *fcnm = "cmopad_uniqueEuler\0";

    *alpha = numpy_mod(*alpha, 2.0*M_PI);

    if (0.5*M_PI < *alpha && *alpha <= M_PI)
    {
        *alpha = M_PI - *alpha;
        *beta  = *beta + M_PI;
        *gamma = 2.0*M_PI - *gamma;
    }
    else if(M_PI < *alpha && *alpha <= 1.5*M_PI)
    {
        *alpha = *alpha - M_PI;
        *gamma = M_PI - *gamma;
    }
    else if(1.5*M_PI < *alpha && *alpha <= 2.0*M_PI)
    {
        *alpha = 2.0*M_PI - *alpha;
        *beta  = *beta + M_PI;
        *gamma = M_PI + *gamma;
    }

    *alpha = numpy_mod(*alpha, 2.0*M_PI);
    *beta  = numpy_mod(*beta,  2.0*M_PI);
    *gamma = numpy_mod(*gamma+M_PI, 2.0*M_PI) - M_PI;

    // If dip is exactly 90 degrees, one is still
    // free to choose between looking at the plane from either side.
    // Choose to look at such that beta is in the range [0,180)

    // This should prevent some problems, when dip is close to 90 degrees:
    if (fabs(*alpha - 0.5*M_PI) < 1e-10){*alpha = 0.5*M_PI;}
    if (fabs(*beta - M_PI) < 1e-10){*beta = M_PI;}
    if (fabs(*beta - 2.*M_PI) < 1e-10){*beta = 0.;}
    if (fabs(*beta) < 1e-10){*beta = 0.;}

    if ((fabs(*alpha - 0.5*M_PI) < 1.e-14) && *beta >= M_PI)
    {
        *gamma =-*gamma;
        *beta  = numpy_mod(*beta-M_PI, 2.0*M_PI);
        *gamma = numpy_mod(*gamma+M_PI, 2.0*M_PI) - M_PI;
        if (!(0. <= *beta && *beta < M_PI))
        {
            printf("%s: Warning on bounds on beta!", fcnm);
        }
        if (!(-M_PI <= *gamma && *gamma < M_PI))
        {
            printf("%s: Warning on bounds on gamma!", fcnm);
        }
    }

    if (*alpha < 1e-7)
    {
        *beta = numpy_mod(*beta + *gamma, 2.0*M_PI);
        *gamma = 0.;
    }
    return;
}
//============================================================================//
/*!
 * @brief Ascending argument sort for three numbers.  For more information:
 *        http://stackoverflow.com/questions/4367745/simpler-way-of-sorting-three-numbers
 *
 * @param[in] x       array to sort into ascending[3]
 *
 * @param[out] iperm  permutation such that x3[iperm[:]] is in ascending
 *                    order
 *
 */
static int __cmopad_argsort3(double x[3], int iperm[3])
{
    const char *fcnm = "__cmopad_argsort3\0";
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
//============================================================================//
/*!
 * @brief Determinant of a 3 x 3 matrix
 *
 * @result determinant of a 3 x 3 matrix
 *
 */
static double __cmopad_determinant3x3(double A[3][3])
{
    double determinant;
    determinant = A[0][0]*( A[1][1]*A[2][2] - A[1][2]*A[2][1]) 
                - A[0][1]*( A[1][0]*A[2][2] - A[1][2]*A[2][0])
                + A[0][2]*( A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    return determinant;
}
//============================================================================//
/*!
 * @brief Computes the inverse of a 3 x 3 matrix
 *
 * @param[in,out] Amat   on input the 3 x 3 matrix to invert.
 *                       on output the inverted matrix
 *
 * @result 0 indicates success
 *
 */
static int __cmopad_inv3(double Amat[3][3])
{
    const char *fcnm = "__cmopad_inv3\0";
    double a11, a12, a13, a21, a22, a23, a31, a32, a33, detA, detAi;
    // Determinant of a
    detA = __cmopad_determinant3x3(Amat);
    if (detA == 0.0)
    {
        printf("%s: Error matrix is singular!\n", fcnm);
        return -1;
    }
    detAi = 1.0/detA;
    a11 = Amat[0][0]; a12 = Amat[0][1]; a13 = Amat[0][2];
    a21 = Amat[1][0]; a22 = Amat[1][1]; a23 = Amat[1][2];
    a31 = Amat[2][0]; a32 = Amat[2][1]; a33 = Amat[2][2];
    // Row 1
    Amat[0][0] = detAi*(a22*a33 - a32*a23);
    Amat[0][1] = detAi*(a13*a32 - a12*a33);
    Amat[0][2] = detAi*(a12*a23 - a13*a22);
    // Row 2
    Amat[1][0] = detAi*(a23*a31 - a21*a33);
    Amat[1][1] = detAi*(a11*a33 - a13*a31);
    Amat[1][2] = detAi*(a13*a21 - a11*a23);
    // Row 3
    Amat[2][0] = detAi*(a21*a32 - a22*a31);
    Amat[2][1] = detAi*(a12*a31 - a11*a32);
    Amat[2][2] = detAi*(a11*a22 - a12*a21); 
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the trace (sum of diagonal) of a 3x3 matrix
 *
 * @param[in] M     matrix of which to compute trace
 *
 * @result trace of matrix
 *
 */
static double __cmopad_trace3(double M[3][3])
{
    double trace;
    trace = M[0][0] + M[1][1] + M[2][2]; 
    return trace;
}
//===========================================================================//
/*!  
 * @brief Calculates all eigenvalues and eigenvectors (job = 0) of a
 *        symmetric 3 x 3 matrix using LAPACK. 
 *
 * @param[in,out] a  on input matrix of which to calculate eigenvalues
 *                   and optionally eigenvectors.
 *                   on output if the eigenvectors are desired they are
 *                   stored on a 
 *
 * @param[in] job    = 0 calculate eigenvectors, otherwise, just eigenvalues
 * 
 * @param[out] eigs  eigenvalues of matrix a [3]
 *
 * @result 0 indicates success
 * 
 * @author B. Baker (ISTI)
 * @date April 2016
 */
static int __cmopad_Eigs3x3(double a[3][3], int job, double eigs[3])
{
    const char *fcnm = "__cmopad_Eigs3x3\0";
    double avec[9];
    int n = 3;
    int i, j, indx, ierr;
    int lda = n;
    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            // fill lower half
            if (i != j)
            {
                indx = i * lda + j;
                avec[indx] = a[j][i];
            }
            indx = j * lda + i;
            avec[indx] = a[i][j];
        }
    }
    // Calculate eigenvalues/eigenvectors
    if (job != 0)
    {
        ierr = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'L', n, avec, lda, eigs);
    }
    else
    {
        ierr = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', n, avec, lda, eigs);
    }
    if (ierr != 0)
    {
        printf("%s: Error computing eigenvalues\n", fcnm);
        ierr = 1;
    }
    // Copy eigenvectors back onto A
    if (job == 0 && ierr == 0)
    {
        indx = 0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                indx = j * lda + i;
                a[i][j] = avec[indx];
            }
        }
    }
    return 0;
}
/*!
 * @brief Function for swapping to two numbers a, b
 *
 * @param[in,out] a    on input a = a, on output a = b
 * @param[in,out] b    on input b = b, on output b = a
 *
 */
static void __cmopad_swap8(double *a, double *b)
{
    double temp;
    temp = *a;
    *a = *b;
    *b = temp;
    return;
}
