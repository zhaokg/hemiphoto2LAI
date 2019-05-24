#pragma once
#include  "abc_mem.h"
#include  "abc_common.h"
#include  "mex.h"
#include  "hemiphoto.h"
extern bool  get_options(OPT * _restrict opt,int32_t nrhs,const mxArray * _restrict prhs[],MemPointers * _restrict MEM);
extern mxArray * allocate_output(RES *  _restrict result,int32_t N,int32_t * nptr);
