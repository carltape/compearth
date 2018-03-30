#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"
static void swap(double *a, double *b);
static void sort3descend(const double *__restrict__ lam3,
                         double *__restrict__ lam3sort);
/*!
 * @brief Sorts eigenvalues as lam1 >= lam2 >= lam3.
 *
 * @param[in] n         Number of eigenvalue triplets.
 * @param[in] lam       The eigenvalues to sort.  This is an array
 *                      of dimension [3 x n] with leading dimension 3.
 *
 * @param[out] lamSort  Eigenvalue triplets sorted such that each
 *                      triplet is in descending order.  This is
 *                      an array of dimension [3 x n] with leading
 *                      dimension 3.
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 * 
 */
int compearth_lamsort(const int n, const double *__restrict__ lam,
                      double *__restrict__ lamSort)
{
    int i;
    if (n < 1 || lam == NULL || lamSort == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: No eigentriples\n", __func__);}
        if (lam == NULL){fprintf(stderr, "%s: lam is NULL\n", __func__);}
        if (lamSort == NULL)
        {
            fprintf(stderr, "%s: lamSort is NULL\n", __func__);
        }
        return -1;
    }
/*
    if (n == 3)
    {
        fprintf(stderr, "%s: lam is 3x3 make sure each column is a lambda vector\n",
               __func__);
    }
*/
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
