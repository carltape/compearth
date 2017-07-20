#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

int check_Udetcheck(void);

int main( )
{
    int ierr;
    ierr = check_Udetcheck();
    if (ierr != 0){printf("failed Udetcheck\n"); return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}

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
        }
    } 
    return EXIT_SUCCESS;
}
