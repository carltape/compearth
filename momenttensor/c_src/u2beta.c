#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes beta lune colatitude from u.  This applies the inverse
 *        of Equation 24a of Tape and Tape 2015 
 *
 * @param[in] n          number of points
 * @param[in] maxit      max number of iterations
 * @param[in] linvType   if =1 -> Newton-Rhapson iteration.
 *                       Otherwise -> Halley iteration (default).
 * @param[in] u          u coordinates \f$ [0, 3*pi/4] \f$ [n]
 * @param[in] tol        convergence achieved when:
 *                         fabs(u - 0.75*beta - 0.5*sin2b + 0.0625*sin4b) < tol
 *
 * @param[out] beta      colalitudes \f$ [0, \pi] \f$ [n]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under MIT
 *
 * @bug there is an issue with small angle approximations near the endpoints
 *
 */
int compearth_u2beta(const int n,
                     const int maxit,
                     const int linvType,
                     const double *__restrict__ u,
                     const double tol,
                     double *__restrict__ beta)
{
    const char *fcnm = "compearth_u2beta\0";
    const int nb = 16;
    const double betas[16] = {0.0000000000, 0.1682996064, 0.3365992129,
                              0.5048988193, 0.6731984258, 0.8414980322,
                              1.0097976387, 1.1780972451, 1.3463968515,
                              1.5146964580, 1.6829960644, 1.8512956709,
                              2.0195952773, 2.1878948838, 2.3561944902,
                              2.3561944902};
    const double us[16] = {0.0000000000, 0.0000532865,
                           0.0016375032, 0.0116225691,
                           0.0445525969, 0.1203598608,
                           0.2579993274, 0.4675195432,
                           0.7439913014, 1.0661325471,
                           1.4006252490, 1.7107983457,
                           1.9665451937, 2.1518309406,
                           2.2671458676, 100.0};
    double absBetaNew, betaNew, betaWork, c, cos2b, cos4b,
           den, dfdx, dfdx2, f, sin2b, sin4b;
    int i, ierr, indx, k;
    bool lconv;
    const double lowerBound = 0.0 + 0.144;
    const double upperBound = 0.75*M_PI - 0.144;
    ierr = 0;
    betaWork = 0.0;
#ifdef USE_OPENMP
  #pragma omp parallel for \
  firstprivate(betaWork), \
  private (betaNew, c, cos2b, cos4b, den, dfdx, dfdx2,  \
           f, i, indx, k, lconv, sin2b, sin4b)  \
  shared (beta, betas, fcnm, lowerBound, u, upperBound ) \
  reduction (max:ierr) \
  default (none)
#endif
    for (i=0; i<n; i++)
    {
        // Edge cases
        ierr = 0;
        if (u[i] < lowerBound || fabs(u[i] - lowerBound) < tol/2.0)
        {
            beta[i] = fmax(0, u[i]); //TODO need small angle approximation
//compearth_beta2u(1, &beta[i], &f);
//printf("diff: %e %e\n", f - u[i], u[i]);
 //           continue;
        }
        if (u[i] > upperBound || fabs(u[i] - upperBound) < tol/2.0)
        {
            // TODO need small angle approximation
            beta[i] = fmax(M_PI, M_PI - fabs(us[nb-2] - u[i]));
//compearth_beta2u(1, &beta[i], &f);
//printf("diff: %e %f\n", f - u[i], u[i]);
//            continue;
        }
        // Bracket the solution 
        indx = 0;
        for (k=0; k<nb-1; k++)
        {
            if (us[k] < u[i] && u[i] < us[k+1])
            {
                indx = k;
                break;
            }
        }
        if (k < nb)
        {
            // try to take closer interval
            if (k < nb - 1)
            {
                if (fabs(us[k+1] - u[i]) < fabs(us[k] - u[i])){indx = k + 1;}
            }
            betaWork = betas[indx];
        }
        else
        {
            printf("%s: Failed to bracket solution %d\n", fcnm, i);
            betaWork = 0.5*(upperBound + lowerBound);
        }
        // Begin Newton-Rhapson iteration
        lconv = false;
        if (linvType == 1)
        {
            for (k=0; k<maxit; k++)
            {
                //f = tape2015_beta2u(betaWork) - u;
                compearth_beta2u(1, &betaWork, &f);
                f = f - u[i];
                // Primary convergence check
                if (fabs(f) < tol)
                {
                    lconv = true;
                    beta[i] = betaWork;
                    break;
                }
                cos2b = cos(2.0*betaWork);
                cos4b = cos(4.0*betaWork);
                dfdx = 0.75 - cos2b + 0.25*cos4b;
                // Avoid a breakdown in Newton Rhapson
                if (fabs(dfdx) < 1.e-15)
                {
                    printf("%s: Division by zero\n", fcnm);
                    lconv = false;
                    break;
                }
                betaNew = betaWork - f/dfdx; // Newton-Rhapson iteration
                absBetaNew = fabs(betaNew);
                if (fabs(betaNew - betaWork)/absBetaNew < tol)
                {
                    lconv = true;
                    beta[i] = betaWork;
                    break; 
                }
                betaWork = betaNew;
            } // Loop on iterations
        }
        // Begin Halley iteration
        else
        {
            for (k=0; k<maxit; k++)
            {
                sin2b = sin(2.0*betaWork);
                sin4b = sin(4.0*betaWork);
                f = 0.75*betaWork - 0.5*sin2b + 0.0625*sin4b - u[i];
                // Primary convergence check
                if (fabs(f) < tol)
                {
                    lconv = true;
                    beta[i] = betaWork;
                    break;
                }
                cos2b = cos(2.0*betaWork);
                cos4b = cos(4.0*betaWork);
                dfdx = 0.75 - cos2b + 0.25*cos4b;
                dfdx2 = 2.0*sin2b - sin4b;
                den = 2.0*dfdx*dfdx - f*dfdx2;
                if (fabs(den) < 1.e-15)
                {
                    printf("%s: Division by zero\n", fcnm);
                    lconv = false;
                    break;
                }
                c = 2.0*f*dfdx/den;
                if (fabs(dfdx) < tol)
                {
                   lconv = true;
                   beta[i] = betaWork;
                   break;
                }
                betaNew = betaWork - c; // Halley iteration
                // Require iteration making good progress
                absBetaNew = fabs(betaNew);
                if (fabs(betaNew - betaWork)/absBetaNew < tol)
                {
                    lconv = true;
                    beta[i] = betaWork;
                    break;
                }
                betaWork = betaNew;
            }
        } // End check on halley vs newton
        // check convergence
        if (!lconv)
        {
            compearth_beta2u(1, &betaWork, &f);
            printf("%s: Failure to converge after %d iterations %d %f %f\n",
                   fcnm, maxit, i, f, u[i]);
            ierr = 1;
            beta[i] = 0.0;
        }
        //compearth_beta2u(1, &beta[i], &f);
        //printf("solved: diff: %e %e\n", f - u[i], u[i]);
    }
    return ierr;
}
