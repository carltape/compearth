#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"
static void swap(double *a, double *b);
static void sort3descend(const double *__restrict__ lam3,
                         double *__restrict__ lam3sort);
/*!
 * @brief Sorts eigenvalues as lam1 >= lam2 >= lam3
 *
 * @param[in] n         number of eigenvalue triplets
 * @param[in] lam       the eigenvalues to sort [3*n]
 *
 * @param[out] lamSort  eigenvalue triplets sorted such that each
 *                      triplet is in descending order [3*n]
 *
 * @result 0 indicates success
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 * 
 */
int compearth_lamsort(const int n, const double *__restrict__ lam,
                      double *__restrict__ lamSort)
{
    const char *fcnm = "compearth_lamsort\0";
    int i;
    if (n < 1 || lam == NULL || lamSort == NULL)
    {
        if (n < 1){printf("%s: No eigentriples\n", fcnm);}
        if (lam == NULL){printf("%s: lam is NULL\n", fcnm);}
        if (lamSort == NULL){printf("%s: lamSort is NULL\n", fcnm);}
        return -1;
    }
    if (n == 3)
    {
        printf("%s: lam is 3x3 make sure each column is a lambda vector\n",
               fcnm);
    }
    for (i=0; i<n; i++)
    {
        sort3descend(&lam[3*i], &lamSort[3*i]);
    }
    return 0;
}

static void sort3descend(const double *__restrict__ lam3,
                         double *__restrict__ lam3sort)
{
    double a, b, c;
    a = lam3[0];
    b = lam3[1];
    c = lam3[2];
    if (a > c){swap(&a, &c);}
    if (a > b){swap(&a, &b);}
    // now a is the biggest element - order b and c
    if (b > c){swap(&b, &c);}
    // put into descending order
    lam3sort[0] = c;
    lam3sort[1] = b;
    lam3sort[2] = a;
    return;
}

static void swap(double *a, double *b)
{
    double temp;
    temp = *a;
    *a = *b;
    *b = temp;
    return;
}
