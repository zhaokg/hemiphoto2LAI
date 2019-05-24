#pragma once
#include "abc_000_macro.h"
#include "float.h"
#define FLOAT_TYPE 8   
#if FLOAT_TYPE==4
#define FLOAT float
#define FLOAT_MAX FLT_MAX
#define FLOAT_MIN -FLT_MAX
#define FLOAT_EPSILON FLT_EPSILON
#elif FLOAT_TYPE==8
#define FLOAT double
#define FLOAT_MAX DBL_MAX
#define FLOAT_MIN -DBL_MAX
#define FLOAT_EPSILON DBL_EPSILON
#endif
#define PI   3.141592653589793
#define ISNAN(f) (f) !=(f)
#define ISINF(f) (f) >FLOAT_MAX||(f)< FLOAT_MIN
#include <inttypes.h>
#define rFLOAT		register FLOAT
#define rFLT		register float
#define rDBL		register double
#define rINT		register int32_t
#define rUINT		register uint32_t
#define rINT16		register int16_t
#define rUINT16		register uint16_t
#define rUINT8		register uint8_t
#define rINT8		register int8_t
typedef FLOAT     * _restrict FLOATPTR;
typedef float   *   _restrict FLTPTR;
typedef double   *  _restrict DBLPTR;
typedef int32_t *   _restrict INTPTR;
typedef uint32_t*   _restrict UINTPTR;
typedef int16_t  *  _restrict INT16PTR;
typedef uint16_t  * _restrict UINT16PTR;
typedef int8_t  *   _restrict INT8PTR;
typedef uint8_t *   _restrict UINT8PTR;
#define rFLOATPTR	register  FLOAT * _restrict
#define rDBLPTR		register  double * _restrict
#define rFLTPTR		register  float * _restrict
#define rINT8PTR    register  int8_t *_restrict
#define rUINT8PTR   register  uint8_t *_restrict
#define rINT16PTR   register  int16_t *_restrict
#define rUINT16PTR  register  uint16_t *_restrict
#define rINTPTR		register  int32_t * _restrict
#define rUINTPTR	register  uint32_t * _restrict
