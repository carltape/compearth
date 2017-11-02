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

#endif /*  COMPEARTH_CONSTANTS_H__ */
