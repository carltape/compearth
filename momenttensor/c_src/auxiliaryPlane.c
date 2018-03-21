#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

#pragma omp declare simd
static void strikeDip(const double nRad, const double eRad, const double uRad,
                      double *strikeDeg, double *dipDeg);

/*!
 * @brief Computes the strike, dip, and rake of the auxiliary nodal plane.
 *
 * @param[in] nmt    Number of moment tensors.
 * @param[in] s1     Strike angles (degrees) for plane 1.  These are measured
 *                   positive east from north and should be in the range of 
 *                   \f$ \phi \in [0, 360] \f$.  This is an array of dimension
 *                   [nmt].
 * @param[in] d1     Dip angles (degrees) for plane 1.  These are measured
 *                   positive down from the horizontal and should be in the
 *                   range \f$ \delta \in [0, 90 \f$.  This is an array
 *                   of dimension [nmt].
 * @param[in] r1     Rake angles (degrees) for plane 1.  This should be in 
 *                   range \f$ \lambda \in [-180, 180] \f$.  This is an array
 *                   of dimension [nmt].
 * @param[out] s2    Strike angles (degrees) for plane 2.  These are measured
 *                   positive east from north and should be in the range of
 *                   \f$ \phi \in [0, 360] \f$.  This is an array of
 *                   dimension [nmt].
 * @param[out] d2    Dip angles (degrees) for plane 2.  These are measured
 *                   positive down from the horizontal and should be in the
 *                   range \f$ \delta \in [0, 90] \f$.  This is an array
 *                   of dimension [nmt].
 * @param[out] r2    Rake angles (degrees) for plane 2.  This should be in 
 *                   range \f$ \lambda \in [-180, 180] \f$.  This is an array
 *                   of dimension [nmt].
 *
 * @result 0 indicates success.
 *
 * @author Andy Michaels and Oliver Boyd updated by Ben Baker
 *
 * @copyright MIT
 *
 * @date August 2017
 *
 */
int compearth_auxiliaryPlane(const int nmt,
                             const double *__restrict__ s1,
                             const double *__restrict__ d1,
                             const double *__restrict__ r1,
                             double *__restrict__ s2,
                             double *__restrict__ d2,
                             double *__restrict__ r2)
{
    const double toRad = M_PI/180.0;
    const double toDeg = 180.0/M_PI;
    double cosz, cosz2, cosz3, hNorm, h1, h2, h3, n1, n2, n3,
           sinz, sinz2, sinz3, sl1, sl2, sl3, z, z2, z3;
    int i;
    if (nmt < 1 ||
        s1 == NULL || d1 == NULL || r1 == NULL ||
        s2 == NULL || d2 == NULL || r2 == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: No points\n", __func__);}
        if (s1 == NULL){fprintf(stderr, "%s: s1 is NULL\n", __func__);}
        if (d1 == NULL){fprintf(stderr, "%s: d1 is NULL\n", __func__);}
        if (r1 == NULL){fprintf(stderr, "%s: r1 is NULL\n", __func__);}
        if (s2 == NULL){fprintf(stderr, "%s: s2 is NULL\n", __func__);}
        if (d2 == NULL){fprintf(stderr, "%s: d2 is NULL\n", __func__);}
        if (r2 == NULL){fprintf(stderr, "%s: r2 is NULL\n", __func__);}
        return -1;
    }
    for (i=0; i<nmt; i++)
    {
        // Advance the strike 90 degrees
        z  = (s1[i] + 90.0)*toRad; // strike +90 degrees
        z2 = d1[i]*toRad; // dip
        z3 = r1[i]*toRad; // rake
        // Compute the slip vector in plane 1
        cosz  = cos(z);
        cosz2 = cos(z2);
        cosz3 = cos(z3); 
        sinz  = sin(z);
        sinz2 = sin(z2);
        sinz3 = sin(z3);
        // Jost and Herrmann Eqn 15 except sl1 and sl3 are negated. 
        sl1 =-cosz3*cosz - sinz3*sinz*cosz2;
        sl2 = cosz3*sinz - sinz3*cosz*cosz2;
        sl3 = sinz3*sinz2; 
        // Compute normal vector to plane 1.  Jost and Herrmann Eqn 16
        // except n1 and n3 are negated.
        n1 = sinz*sinz2;
        n2 = cosz*sinz2;
        n3 = cosz2;
        // Compute slip vector of plane 2.  This looks like:
        //   h = (1,0,0) x (sl1, sl2, sl3)
        // I think this produces the directions of strike and slip.
        hNorm = 1.0/hypot(sl2, sl1); // Letting h3 = 0.0
        h1 =-sl2*hNorm;
        h2 = sl1*hNorm;
        h3 = 0.0;
        // Dot strike and slip into the fault plane (normal to fault plane 1)
        z = h1*n1 + h2*n2 + h3*n3;
        // z/hypot(h1, h2);
        z = fmin(1.0, fmax(-1.0, z)); // Strictly enforce bounds
        // Compute the angle between the measured fault plane and directions
        // of slip and dip - Aki and Richards pg 101 Figure 4.13
        z = acos(z);
        // Compute the strike and dip angles from the slip
        strikeDip(sl2, sl1, sl3, &s2[i], &d2[i]);
        r2[i] = z*toDeg;
        if (sl3 <= 0.0){r2[i] =-r2[i];}
    }
    return 0;
}
/*!
 * @brief Finds strike and dip of plane given normal vector having components
 *        n, e, and u. 
 *
 * @author Andy Michaels and Oliver Boyd updated by Ben Baker
 *
 * @copyright MIT
 *
 * @date August 2017
 *
 */
static void strikeDip(const double nRad, const double eRad, const double uRad,
                      double *strikeDeg, double *dipDeg)
{
    const double toDeg = 180./M_PI; 
    double e, u, n, strike, x;
    n = nRad;
    e = eRad;
    u = uRad;
    if (u < 0.0)
    {
        n =-n;
        e =-e;
        u =-u;
    }
    strike = atan2(e, n)*toDeg; // This is in range [-180,180]
    strike = strike - 90.0; // This is in range [-270,90]
    // This loop is pointless
    //while (strike >= 360.0)
    //{
    //    strike = strike - 360.0;
    //}
    // This loop could execute at most once
    //while strike < 0
    //{
    //    strike = strike + 360.0;
    //}
    if (strike < 0.0){strike = strike + 360.0;}
    *strikeDeg = strike;
    x = hypot(n, e); //x = sqrt(n*n + e*e);
    *dipDeg = atan2(x, u)*toDeg;
    return;
}
/*
function [strike, dip, rake] = AuxPlane(s1,d1,r1);
%function [strike, dip, rake] = AuxPlane(s1,d1,r1);
% Get Strike and dip of second plane, adapted from Andy Michael bothplanes.c
r2d = 180/pi;

z = (s1+90)/r2d;
z2 = d1/r2d;
z3 = r1/r2d;
% slick vector in plane 1 
sl1 = -cos(z3).*cos(z)-sin(z3).*sin(z).*cos(z2);
sl2 = cos(z3).*sin(z)-sin(z3).*cos(z).*cos(z2);
sl3 = sin(z3).*sin(z2);
[strike, dip] = strikedip(sl2,sl1,sl3);

n1 = sin(z).*sin(z2);  % normal vector to plane 1 
n2 = cos(z).*sin(z2);
n3 = cos(z2);
h1 = -sl2; % strike vector of plane 2 
h2 = sl1;
% note h3=0 always so we leave it out 

z = h1.*n1 + h2.*n2;
z = z./sqrt(h1.*h1 + h2.*h2);
z = acos(z);

rake = zeros(size(strike));
j = find(sl3 > 0);
rake(j) = z(j)*r2d;
j = find(sl3 <= 0);
rake(j) = -z(j)*r2d;


function [strike, dip] = strikedip(n, e, u)
%function [strike, dip] = strikedip(n, e, u)
%       Finds strike and dip of plane given normal vector having components n, e, and u
%

% Adapted from Andy Michaels stridip.c

r2d = 180/pi;

j = find(u < 0);
n(j) = -n(j);
e(j) = -e(j);
u(j) = -u(j);

strike = atan2(e,n)*r2d;
strike = strike - 90;
while strike >= 360
        strike = strike - 360;
end
while strike < 0
        strike = strike + 360;
end

x = sqrt(n.^2 + e.^2);
dip = atan2(x,u)*r2d;

*/
