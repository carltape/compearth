%module compearth
%{
#include "compearth.h"
extern int compearth_CMT2omega(const int nmt1, const double * M1,
                        const int nmt2, const double * M2,
                        double * omega);
extern int compearth_standardDecomposition(const int nmt,
                                    const double * M,
                                    enum compearthCoordSystem_enum basis,
                                    double * M0,
                                    double * Mw,
                                    double * fp1,
                                    double * fp2,
                                    double * pAxis,
                                    double * bAxis,
                                    double * tAxis,
                                    double * isoPct,
                                    double * devPct,
                                    double * dcPct,
                                    double * clvdPct);
%}

extern int compearth_CMT2omega(const int nmt1, const double * M1,
                        const int nmt2, const double * M2,
                        double * omega);
extern int compearth_standardDecomposition(const int nmt,
                                    const double * M,
                                    enum compearthCoordSystem_enum basis,
                                    double * M0,
                                    double * Mw,
                                    double * fp1,
                                    double * fp2,
                                    double * pAxis,
                                    double * bAxis,
                                    double * tAxis,
                                    double * isoPct,
                                    double * devPct,
                                    double * dcPct,
                                    double * clvdPct);
