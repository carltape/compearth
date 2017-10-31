#ifndef COMPEARTH_CONSTANTS_H__
#define COMPEARTH_CONSTANTS_H__ 1
#include <stdbool.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_SQRT2
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif
#ifndef M_SQRT1_2
# define M_SQRT1_2      0.70710678118654752440 /* 1/sqrt(2) */
#endif
#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

enum compearthCoordSystem_enum
{
    CE_UNKNOWN_COORDS = 0, /*!< Unkown coordinate system. */
    CE_USE = 1,    /*!< Up, south, east (G-CMT) */
    CE_NED = 2,    /*!< North, east, down (Aki and Richards, 1980 pg 118) */
    CE_NWU = 3,    /*!< North, west, up */
    CE_ENU = 4,    /*!< East, north, up */
    CE_SEU = 5     /*!< South, east, up */
};

enum magType_enum
{
    CE_UNKNOWN_MW = 0,    /*!< Unknown Mw scale. */
    CE_KANAMORI_1978 = 1, /*!< Mw frrom Kanamori Mw = (2/3)*log10(M0) + k; */
    CE_HARVARD_CMT = 2    /*!< Mw from Harvard CMT: (2/3)*(log10(M0) - 16.1) */
};


#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif
enum ceNormType_enum
{
    CE_UNKNOWN_NORM = 0,          /*!< Unknown norm. */
    CE_TWO_NORM = 2,              /*!< \$ L_2 = \sqrt{\sum_i x_i^2} \$ norm */
    CE_ONE_NORM = 1,              /*!< \$ L_1 = \sum_i |x_i| \$ norm */
    CE_P_NORM = 3,                /*!< \$ L_p 
                                    = \left (
                                        \sum_i |x_i|^p \right )^{1/p}
                                      \right ) \$ norm */
    CE_INFINITY_NORM = 4,          /*!< \$ L_\infty = max |x| \$ */
    CE_NEGATIVE_INFINITY_NORM = 5  /*!< \$ L_{-\infty} = min |x| \$ */
};
#ifdef __clang__
#pragma clang diagnostic pop
#endif

enum ceOrthoType_enum
{
    CE_NO_ORTH = 0,   /*!< Don't perform a reorthgonalization. */
    CE_ORTH_SVD = 1,       /*!< Orthgonalizes with SVD. */
    CE_ORTH_TAPE2012 = 2,  /*!< Orthgonalizes with Tape and Tape 2012c Appendix E. */
    CE_ORTH_QUAT = 3       /*!< Orthoganlizes with quaternions. */
};

#endif /*  COMPEARTH_CONSTANTS_H__ */
