#include "gfun.h"
#include "abc_common.h"
#include "abc_gaussquad.h"
FLOATPTR lroots;
FLOATPTR weight;
int32_t  NQ;
int32_t strcicmp(char const * _restrict a,char const * _restrict b);
int32_t gfuncSelect( char *name,char **nameList,int32_t N)
{
	for (rINT i=0; i < N; i++)
	{
		if (strcicmp(name,nameList[i])==0)
			return i;
	}
	return -5;
}
void selectGfunc(char *name,PAR_STRUCT * par)
{
	if (strcicmp(name,"uni")==0)
	{
		par->gfun=gfunc_uni;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"sph")==0)
	{
		par->gfun=gfunc_sph;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"hor")==0)
	{
		par->gfun=gfunc_hor;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"vtc")==0)
	{
		par->gfun=gfunc_vtc;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"erc")==0)
	{
		par->gfun=gfunc_erc;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"pln")==0)
	{
		par->gfun=gfunc_pln;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"plg")==0)
	{
		par->gfun=gfunc_plg;
		par->Npar=0; 	
	}
	else if (strcicmp(name,"ext")==0)
	{
		par->gfun=gfunc_ext;
		par->Npar=0; 
	}
	else if (strcicmp(name,"bet")==0)
	{
		par->gfun=gfunc_bet;
		par->Npar=2;
		par->Gparam_init[0]=par->Gparam_init[1]=2;
		par->boundG[0]=1; par->boundG[1]=60;
		par->boundG[2]=1; par->boundG[3]=60;
	}
	else if (strcicmp(name,"elt")==0)
	{
		par->gfun=gfunc_elt;
		par->Npar=2;
		par->Gparam_init[0]=par->Gparam_init[1]=.5;
		par->boundG[0]=0; par->boundG[1]=1.;
		par->boundG[2]=0; par->boundG[3]=PI/2;
	}
	else if (strcicmp(name,"rgm")==0)
	{
		par->gfun=gfunc_rgm;
		par->Npar=1;
		par->Gparam_init[0]=.1;
		par->boundG[0]=-0.40; par->boundG[1]=0.6;
	}
	else if (strcicmp(name,"dic")==0)
	{
		par->gfun=gfunc_dic;
		par->Npar=1;
		par->Gparam_init[0]=.1;
		par->boundG[0]=-1; par->boundG[1]=1;
	}
	else if (strcicmp(name,"lin")==0)
	{
		par->gfun=gfunc_lin;
		par->Npar=1;
		par->Gparam_init[0]=.1;
		par->boundG[0]=0; par->boundG[1]=1;
	}
	else if (strcicmp(name,"es1")==0)
	{
		par->gfun=gfunc_es1;
		par->Npar=1;
		par->Gparam_init[0]=1;
		par->boundG[0]=0; par->boundG[1]=50;
	}
	else if (strcicmp(name,"es2")==0)
	{
		par->gfun=gfunc_es2;
		par->Npar=1;
		par->Gparam_init[0]=1;
		par->boundG[0]=0; par->boundG[1]=50;
	}
	else if (strcicmp(name,"res")==0)
	{
		par->gfun=gfunc_res;
		par->Npar=1;
		par->Gparam_init[0]=1;
		par->boundG[0]=0; par->boundG[1]=50;
	}
	else if (strcicmp(name,"sts")==0)
	{
		par->gfun=gfunc_sts;
		par->Npar=1;
		par->Gparam_init[0]=.1;
		par->boundG[0]=0; par->boundG[1]=1;
	}
	else if (strcicmp(name,"lan")==0)
	{
		par->gfun=gfunc_lan;
		par->Npar=1;
		par->Gparam_init[0]=.1;
		par->boundG[0]=0; par->boundG[1]=1;
	}	 
	else if (strcicmp(name,"vhf")==0)
	{
		par->gfun=gfunc_vhf;
		par->Npar=2;
		par->Gparam_init[0]=.2;
		par->Gparam_init[1]=.3;
		par->boundG[0]=-1; par->boundG[1]=1;
		par->boundG[2]=-1; par->boundG[3]=1;
	}
}
void gfunc_elt(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT e=gParam[0];
	FLOAT m=gParam[1];
	FLOAT D;
	{
		rFLOAT n=asin(e*cos(m));
		rFLOAT v=asin(e*sin(m));
		D=e/(sin(m)*log((cos(n)+sin(v))/(cos(v) - sin(n)))+cos(m)*(n - v));
	}
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=e*cos(T - m);
			g=1 - g*g;
			g=sin(T)*D/sqrt(g);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=(e*cos(T - m));
			g=1 - g*g;
			g=sin(T)*D/sqrt(g);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3);
	}
}
void gfunc_bet(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT u=gParam[0];
	FLOAT v=gParam[1];
	FLOAT D;
	D=2/PI *tgamma(u+v)/tgamma(u)/tgamma(v);
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ;j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=2 * T/PI;
			g=pow(1 - g,u - 1)*pow(g,v - 1);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=2 * T/PI;
			g=pow(1 - g,u - 1)*pow(g,v - 1);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
void gfunc_es1(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT x=gParam[0];
	rFLOAT x2=x*x;
	rFLOAT f=1/(x+1.702*pow(x+1.12,-0.708));
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		rFLOAT c=cos(t);
		rFLOAT s=sin(t);
		G[i]=sqrt(x2*c*c+s*s)*f;
	}
}
void gfunc_es2(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT x=gParam[0];
	rFLOAT x2=x*x;
	rFLOAT f=1/(x+1.774*pow(x+1.182,-0.733));
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		rFLOAT c=cos(t);
		rFLOAT s=sin(t);
		G[i]=sqrt(x2*c*c+s*s)*f;
	}
}
void gfunc_res(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT x=gParam[0]; 
	FLOAT A;
	if (fabs(x - 1) <  1e-10)
		A=2;	
	else if (x < 1)
	{
		rFLOAT e;
		e=sqrt(1 - x*x);
		A=x+asin(e)/e;
	}
	else if (x>1)
	{
		rFLOAT e;
		e=sqrt(1 - 1/(x*x));
		A=x+log((1+e)/(1-e))/( 2*e*x);
	}
	FLOAT D;
	D=2*x*x*x/A;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			rFLOAT cosTheta=sin(T);
			rFLOAT sinTheta=cos(T);
			g=(sinTheta*sinTheta+x*x*cosTheta*cosTheta);
			g=cosTheta/(g*g);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			rFLOAT cosTheta=sin(T);
			rFLOAT sinTheta=cos(T);
			g=(sinTheta*sinTheta+x*x*cosTheta*cosTheta);
			g=cosTheta/(g*g);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
#define NBIN_GUNFC 100
void gfunc_vhf(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT a=gParam[0];
	FLOAT b=gParam[1];
	{
		rFLOAT sum=fabs(a)+fabs(b);
		if ((sum - 1) > -0.15)
			a=a/(sum+0.15),
			b=b/(sum+0.15);
	}
	FLOAT gDATA[NBIN_GUNFC+2];
	int32_t  Nbin_gfun=NBIN_GUNFC;
	rFLOAT   bin=(PI/2)/Nbin_gfun;
	{
		for (rINT i=0; i <(Nbin_gfun+1); i++)
		{
			rFLOAT t=i*bin;
			rFLOAT x,y;
			rFLOAT c,s,c2;
			rFLOAT dx;			
			if (t> FLOAT_EPSILON)
			{
				rINT ite=0;	
				do{
					c=cos(x),s=sin(x),c2=cos(x+x);
					dx=(t+t+s*(a+b*c) - x)/(a*c+b*c2 - 1);
					x+=-dx;
					ite++;
					if (ite > 100) break;
				} while ((fabs(dx) > 1e-12));
			}	
			else
			{
				x=0;
				gDATA[i]=(2/(1 - a - b) - 1)/(PI/2);
				continue;
			}
			rFLOAT tmp;
			tmp=2/(1 - a*c - b*c2);
			gDATA[i]=(tmp - 1)/(PI/2);
		}
	}   
	gDATA[Nbin_gfun+1]=gDATA[Nbin_gfun];
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		FLOAT  a=0,b=t0;
		FLOAT  c1=(b - a)*0.5,c2=(b+a)*0.5;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			rFLOAT  f_idx=T/bin;
			int32_t i_idx=f_idx;
			rFLOAT  frac=f_idx - i_idx;
			g=(1 - frac)*gDATA[i_idx]+frac*gDATA[i_idx+1];
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			rFLOAT  f_idx=T/bin;
			int32_t i_idx=f_idx;
			rFLOAT  frac=f_idx - i_idx;
			g=(1 - frac)*gDATA[i_idx]+frac*gDATA[i_idx+1];
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3);
	}
}
void gfunc_erc(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=1+cos(2 * T);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=1+cos(2 * T);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
void gfunc_ext(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=1+cos(4 * T);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=1+cos(4 * T);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
void gfunc_hor(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	for (rINT i=0; i < N; i++)
	{
		G[i]=cos(theta[i]);
	}
}
void gfunc_sts(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT x=gParam[0];
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		G[i]=x*cos(t)+(1 - x) * 2/PI*sin(t);
	}
}
void gfunc_lan(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT x=gParam[0];
	rFLOAT x2=(1 - x)*.5;
	for (rINT i=0; i < N; i++)
	{
		G[i]=x*0.5+x2*theta[i];
	}
}
void gfunc_plg(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT e=gParam[0];
	FLOAT m=gParam[1];
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=1 - cos(4 * T);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=1 - cos(4 * T);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
void gfunc_pln(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	FLOAT e=gParam[0];
	FLOAT m=gParam[1];
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=1 - cos(2 * T);
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=1 - cos(2 * T);
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3)*D;
	}
}
void gfunc_rgm(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT XL=gParam[0];
	rFLOAT phi1=0.5 - 0.633*XL - 0.33*XL*XL;
	rFLOAT phi2=(1 - 2 * phi1);
	for (rINT i=0; i < N; i++)
	{
		G[i]=phi1+phi2*cos(theta[i]);
	}
}
void gfunc_dic(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT XL=gParam[0];
	rFLOAT phi1=0.5 - 0.489*XL - 0.11*XL*XL;
	rFLOAT phi2=(1 - 2 * phi1);
	for (rINT i=0; i < N; i++)
	{
		G[i]=phi1+phi2*cos(theta[i]);
	}
}
void gfunc_lin(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT phi1=gParam[0];
	rFLOAT phi2=(1 - 2 * phi1);
	for (rINT i=0; i < N; i++)
	{
		G[i]=phi1+phi2*cos(theta[i]);
	}
}
void gfunc_sph(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	for (rINT i=0; i < N; i++)
	{
		G[i]=.5f;
	}
}
void gfunc_uni(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		FLOAT  cost=cos(t);
		FLOAT  sint=sin(t);
		rFLOAT t0=PI/2 - t;
		double a=0,b=t0;
		double c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT FSUM1=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f,g;
			g=D;
			f=g*cos(T);
			FSUM1+=weight[j] * f;
		}
		FSUM1=c1 * FSUM1*cost;
		a=t0,b=PI/2;
		c1=(b - a)/2,c2=(b+a)/2;
		rFLOAT cott=tan(PI/2 - t);
		rFLOAT FSUM2=0;
		rFLOAT FSUM3=0;
		for (rINT j=0; j < NQ; j++)
		{
			rFLOAT T=c1 * lroots[j]+c2;
			rFLOAT f2,f3,g;
			g=D;
			rFLOAT phi=acos(-cott*tan(PI/2 - T));
			if (phi==phi)
			{
				f2=g*cos(T)*(2/PI*phi - 1);
				f3=g*sin(T)*cos(phi - PI/2);      
			}
			FSUM2+=weight[j] * f2;
			FSUM3+=weight[j] * f3;
		}
		FSUM2=c1*FSUM2*cost;
		FSUM3=c1*FSUM3 * 2/PI*sint;
		G[i]=(FSUM1+FSUM2+FSUM3);
	}
}
void gfunc_vtc(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem)
{
	rFLOAT D=2/PI;
	for (rINT i=0; i < N; i++)
	{
		rFLOAT t=theta[i];
		G[i]=D*sin(t);
	}
}
