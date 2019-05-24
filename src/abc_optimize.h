#pragma once
#include "abc_common.h"
typedef FLOAT(*FUNPTR) (FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern  void  miminize(FUNPTR f,FLOATPTR X,int32_t N,void *par,int32_t length,FLOATPTR mem);
extern  FLOAT rosenbrock(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
