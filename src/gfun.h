#pragma once
#include "hemiphoto.h"
extern FLOATPTR lroots;
extern FLOATPTR weight;
extern int32_t  NQ;
extern int32_t gfuncSelect(char *name,char **nameList,int32_t N);
extern void selectGfunc(char *name,PAR_STRUCT * par);
extern void gfunc_elt(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_bet(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_es1(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_es2(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_res(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_erc(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_ext(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_hor(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_sts(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_lan(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_plg(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_pln(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_rgm(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_dic(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_lin(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_sph(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_uni(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_vtc(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
extern void gfunc_vhf(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
