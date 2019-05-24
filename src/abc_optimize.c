#include "abc_optimize.h"
static INLINE FLOAT dotProdcut(int32_t N,FLOATPTR X,FLOATPTR Y)
{
	rFLOAT asum=0;
	for (rINT i=0; i < N; i++) asum+=X[i] * Y[i];
	return asum;
}
void miminize(FUNPTR f,FLOATPTR X,int32_t N,void *par,int32_t length,FLOATPTR mem)
{
	static const FLOAT INT=0.1f;    
	static const FLOAT EXT=3.0f;    
	static const int32_t MAX=20;    
	static const FLOAT RATIO=10.f;  
	static const FLOAT SIG=0.1f;    
	static const FLOAT RHO=0.05f;   
	static const FLOAT red=1.f;
	FLOAT nan=(FLOAT_MAX*FLOAT_MAX) * 0;
	int32_t i=0L;
	uint8_t ls_failed=0L;
	FLOAT  f0,d0;
	FLOATPTR df0;
	FLOATPTR s; 
	df0=mem;
	s=mem+N;
	f0=f(X,df0,N,par);
	i+=(length<0);            
	d0=0;
	for (rINT i=0; i < N; i++)
	{
		s[i]=-df0[i]; d0+=s[i] * s[i];
	}
	d0=-d0;
	FLOATPTR Xtmp1=mem+2 * N;
	FLOATPTR Xtmp2=mem+3 * N;
	FLOAT x3=red/(1. - d0); 
	FLOATPTR df3=mem+4 * N;
	int32_t LENGTH=abs(length);
	while (i < LENGTH)
	{
		i+=(length>0);
		FLOATPTR X0=X;
		FLOAT    F0=f0;
		FLOATPTR dF0=df0;
		int32_t M;
		M=length > 0 ? MAX : min(MAX,-length - i);
		FLOAT  x2,f2,d2;
		FLOAT  f3,d3; 
		while (1) 
		{
			x2=0; f2=f0; d2=d0;
			char success=0;
			while (!success && M > 0)
			{
				M -=1; i+=(length<0);                         
				for (rINT i=0; i < N; i++) Xtmp1[i]=X[i]+x3*s[i];
				f3=f(Xtmp1,df3,N,par);
				rFLOAT asum=0;
				for (rINT i=0; i < N; i++) asum+=df3[i];
				if (ISNAN(f3)||ISINF(f3)||ISNAN(asum)||ISINF(asum))
				{
					x3=(x2+x3)*0.5;              
					continue;
				}
				success=1;
			}
			if (f3 < F0)
			{
				X0=Xtmp1; F0=f3; dF0=df3;    
			}
			d3=dotProdcut(N,df3,s); 
			if (d3 > SIG*d0||f3 > f0+x3*RHO*d0||M==0) 
			{
				break;
			}
			FLOAT x1=x2,f1=f2,d1=d2;                
			x2=x3; f2=f3; d2=d3;                      
			{
				rFLOAT x2_x1=x2 - x1;
				FLOAT A=6 * (f1 - f2)+3 * (d2+d1)*x2_x1;     
				FLOAT B=3 * (f2 - f1) - (2 * d1+d2)*x2_x1;
				rFLOAT delta=B*B - A*d1*x2_x1;
				if (delta >=0)
					x3=x1 - d1*(x2_x1 *x2_x1)/(B+sqrt(delta)); 
				else
					x3=nan;
			}
			if (ISNAN(x3)||ISINF(x3)||x3 < 0) 
				x3=x2*EXT;                        
			else if (x3 > x2*EXT)                   
				x3=x2*EXT;                        
			else if (x3 < x2+INT*(x2 - x1))      
				x3=x2+INT*(x2 - x1);
		}    
		FLOAT x4,f4,d4;
		x4=x3,f4=f3,d4=d3;
		while ((abs(d3) > -SIG*d0||f3 > f0+x3*RHO*d0) && M > 0) 
		{
			if (d3 > 0||f3 > f0+x3*RHO*d0)                 
				x4=x3,f4=f3,d4=d3;                     
			else
				x2=x3,f2=f3,d2=d3;                     
			if (f4 > f0)
			{
				rFLOAT x4_x2=(x4 - x2);
				x3=x2 - 0.5*d2*(x4_x2*x4_x2)/(f4 -f2 - d2*x4_x2);  
			}
			else
			{
				rFLOAT x4_x2=(x4 - x2);
				rFLOAT A=6 * (f2 - f4)/x4_x2+3 * (d4+d2);   
				rFLOAT B=3 * (f4 - f2) - (2 * d2+d4)*x4_x2;
				rFLOAT delta=B*B - A*d2*(x4_x2*x4_x2);
				if (delta >=0)
				{
					rFLT a=(sqrt(delta) );
					x3=x2+(a-B)/A; 
				}	
				else
					x3=nan;
			}
			if (ISNAN(x3)||ISINF(x3))
			{
				x3=(x2+x4)/2;            
			}
			x3=max(min(x3,x4 - INT*(x4 - x2)),x2+INT*(x4 - x2));  
			for (rINT i=0; i < N; i++) Xtmp2[i]=X[i]+x3*s[i];
			f3=f(Xtmp2,df3,N,par);
			if (f3 < F0)
			{
				X0=Xtmp2; F0=f3; dF0=df3;        
			}
			M=M - 1; i=i+(length<0);             
			d3=dotProdcut(N,df3,s);                        
		}
		if (abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0)
		{
			f0=f3; 
			rFLOAT df3_df03=0;
			rFLOAT df0_sqr_factor=0;
			for (rINT i=0; i < N; i++)
			{
				X[i]+=x3*s[i];   
				df3_df03+=(df3[i] - df0[i]) * df3[i];
				df0_sqr_factor+=df0[i] * df0[i];
			}
			d3=d0;
			d0=0;
			df0_sqr_factor=df3_df03/df0_sqr_factor;
			for (rINT i=0; i < N; i++)
			{
				s[i]=df0_sqr_factor*s[i] - df3[i];
				df0[i]=df3[i];
				d0+=df0[i] * s[i];
			}
			if (d0 > 0) 
			{
				d0=0;
				for (rINT i=0; i < N; i++)
				{
					s[i]=-df0[i];
					d0+=s[i] * s[i];
				}
				d0=-d0;
			}
			x3=x3 * min(RATIO,d3/(d0 - FLOAT_EPSILON));  
			ls_failed=0;                        
		}
		else
		{
			f0=F0;
			d0=0;
			for (rINT i=0; i<N; i++)
			{
				df0[i]=dF0[i];
				X[i]=X0[i];
				s[i]=-X0[i];
				d0+=X0[i] * X0[i];
			}
			d0=-d0;
			if (ls_failed||i > LENGTH)
			{
				break;                       
			}
			x3=1/(1 - d0);
			ls_failed=1;                       
		}
	}
}
FLOAT  rosenbrock(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOAT f=0;
	for (rINT i=0; i < (N-1); i++)
	{
		rFLOAT a=(X[i+1] - X[i] * X[i]);
		rFLOAT b=(1 - X[i]);
		f=f+100 * a*a+b*b;
		df[i]=-400 * X[i] * (X[i+1] - X[i] * X[i]) - 2 * (1 - X[i]);		
	}
	df[N - 1]=0;
	for (rINT i=1; i < N ; i++)
	{		
		df[i]+=200 * (X[i] - X[i - 1] * X[i - 1]);
	}
	return f;
}
