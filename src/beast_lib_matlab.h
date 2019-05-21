#pragma once
#include  "abc_mem.h"
#include  "abc_common.h"
#include  "beast_common.h"
#include  "mex.h"
extern bool  get_options(Options * _restrict opt,const mxArray * _restrict pY,const mxArray * _restrict pOpt,MemPointers * _restrict MEM);
extern mxArray *allocate_output(RESULT *  _restrict matOutput,Options * _restrict opt,int * nptr);
extern bool     get_options_trend(Options * _restrict opt,const mxArray * _restrict pY,const mxArray * _restrict pOpt,MemPointers * _restrict MEM);
extern mxArray *allocate_output_trend(RESULT * _restrict matOutput,Options * _restrict opt,int * nptr);
