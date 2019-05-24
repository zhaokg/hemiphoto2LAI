#include "abc_gaussquad.h"
#include "math.h"
static FLOAT *lroots;
static FLOAT *weight;
static FLOAT *lcoef;
static int32_t N;
static INLINE void lege_init()
{ 
	rINT N2=(N+1)*(N+1);
	for (rINT i=0; i < N2; i++) lcoef[i]=0;
}
static void lege_coef()
{
	rINT N1=N+1;
	*(lcoef+0*N1+0)=1; *(lcoef+1 * N1+1)=1;
	for (rINT n=2; n <=N; n++) {
		*(lcoef+n * N1+0)=-(n - 1) * (*(lcoef+(n-2) * N1+0))/n;
		for (rINT i=1; i <=n; i++)
			*(lcoef+n * N1+i)=((2 * n - 1) *  (*(lcoef+(n-1) * N1+i-1)) 
			- (n - 1) * (*(lcoef+(n - 2) * N1+i )) )/n;
	}
}
static FLOAT lege_eval(int n,FLOAT x)
{
	rINT N1=N+1;
	FLOAT s=*(lcoef+n * N1+n) ;
	for (rINT i=n; i; i--)
		s=s * x+*(lcoef+n * N1+i-1);
	return s;
}
static FLOAT lege_diff(int n,FLOAT x)
{
	return n * (x * lege_eval(n,x) - lege_eval(n - 1,x))/(x * x - 1);
}
void lege_roots(FLOATPTR ROOTS,FLOATPTR WEIGHTS,FLOATPTR COEF_MEM,int32_t knotNumber)
{
	N=knotNumber;
	lroots=ROOTS;
	weight=WEIGHTS;
	lcoef=COEF_MEM;
	lege_init();
	lege_coef();
	FLOAT x,x1;
	for (rINT i=1; i <=N; i++) {
		x=cos(PI * ((N-i+1) - .25)/(N+.5));
		do {
			x1=x;
			x -=lege_eval(N,x)/lege_diff(N,x);
		} while (fdim(x,x1) > 2e-106);
		lroots[i - 1]=x;
		x1=lege_diff(N,x);
		weight[i - 1]=2/((1 - x * x) * x1 * x1);
	}
}
