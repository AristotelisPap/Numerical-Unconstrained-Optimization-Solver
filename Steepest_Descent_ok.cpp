#include<stdio.h>
#include<math.h>

double fun (double x[],int a,int n)
{
    int i;
	double func=0;
	for (i=0;i<n/2;i++)
	{
		func=func+a*pow(x[i+1]-pow(x[i],2),2)+pow(1-x[i],2);	
	}
    return (func);
}

double subgradient_odd (double xi,double xiplus1,int a)
{
    return ((-4)*a*xi*(xiplus1-pow(xi,2))+2*(xi-1));
}

double subgradient_even (double xi,double xiplus1,int a)
{
    return (2*a*(xiplus1-pow(xi,2)));
}

double armijo (double x[],double d[],double gradient[],int a,int n)
{
	int i;
	float theta=0.2,eta=2;
	double lam=1,funpr,fdd=0,funx,fxeta;
	funpr=fun(x,a,n);
	for (i=0;i<n;i++)
	{
		fdd=fdd+gradient[i]*d[i];
		x[i]=x[i]+lam*d[i];
	}
	funx=fun(x,a,n);
	if (fun(x,a,n)<=funpr+theta*lam*fdd)
	{
		for (i=0;i<n;i++)
		{
			x[i]=x[i]+lam*d[i];	
		}
		fxeta=fun(x,a,n);
		if (fxeta<=funpr+eta*theta*lam*fdd)
		{
			lam=2;
			for (i=0;i<n;i++)
			{
				x[i]=x[i]+lam*d[i];	
			}
			fxeta=fun(x,a,n);
		}
		while (fun(x,a,n)<=funpr+eta*theta*lam*fdd)
		{
			lam=eta*lam;
			for (i=0;i<n;i++)
			{
				x[i]=x[i]+lam*d[i];	
			}
			fxeta=fun(x,a,n);
		}
	}
	else
	{
		lam=0.5;
		for (i=0;i<n;i++)
		{
			x[i]=x[i]-lam*d[i];	
		}
		funx=fun(x,a,n);
		while (fun(x,a,n)>funpr+theta*lam*fdd)
		{
			lam=lam/eta;
			for (i=0;i<n;i++)
			{
				x[i]=x[i]-lam*d[i];	
			}
			funx=fun(x,a,n);
		}
	}
	printf("directional=%f,step=%2.18f\n",fdd,lam);
	return (lam);
}



int main()
{
	/* Steepest Descent with Armijo step */
	int a,i,j,choice,n,iter=0;
	printf("Choose the parameters a and n:\n");
	scanf("%d %d",&a,&n);
	double step=1,x[n],xk[n],d[n],gradient[n],fun_eval;
	double epsilon=pow(10,-5),norm=0,ssgrad=0;
	for (i=0;i<n;i++)
	{
		gradient[i]=0;
	}
	for (i=0;i<n;i=i+2)
	{
		xk[i]=-1;
		x[i]=-1;
	}
	for (i=1;i<n;i=i+2)
	{
		xk[i]=-1;
		x[i]=-1;
	}
	fun_eval=fun(x,a,n);
	for (i=0;i<n;i=i+2)
	{
		gradient[i]=subgradient_odd(x[i],x[i+1],a);
		d[i]=-gradient[i];
	}
	for (i=1;i<n;i=i+2)
	{
		gradient[i]=subgradient_even(x[i-1],x[i],a);
		d[i]=-gradient[i];
	}
	for (i=0;i<n;i++)
	{
		printf("gradient%d=%f\n",i+1,gradient[i]);
	}
	for (i=0;i<n;i++)
	{
		ssgrad=ssgrad+pow(gradient[i],2);
	}
	norm=sqrt(ssgrad);
	for (i=0;i<n;i++)
	{
		d[i]=d[i]/norm;
	}
	while (norm>=epsilon)
	{
		step=armijo(x,d,gradient,a,n);
		for (i=0;i<n;i++)
		{
			xk[i]=xk[i]+step*d[i];
			x[i]=xk[i];
			printf("xk%d=%f\n",i+1,xk[i]);
		}
		for (i=0;i<n;i=i+2)
		{
			gradient[i]=subgradient_odd(x[i],x[i+1],a);
			d[i]=-gradient[i];
		}
		for (i=1;i<n;i=i+2)
		{
			gradient[i]=subgradient_even(x[i-1],x[i],a);
			d[i]=-gradient[i];
		}
		ssgrad=0;
		for (i=0;i<n;i++)
		{
			ssgrad=ssgrad+pow(gradient[i],2);
		}
		norm=sqrt(ssgrad);
		for (i=0;i<n;i++)
		{
			d[i]=d[i]/norm;
		}
		iter=iter+1;
		printf("f(x)=%3.8f\n",fun(xk,a,n));
		printf("Iteration %d finished!\n\n",iter);		
	}
	printf("\tAfter %d iter, we get:\n",iter);
	for (i=0;i<n;i++)
	{
		printf("\tx%d,opt=%2.5f\n",i+1,xk[i]);
	}
	printf("\tf(xopt)=%3.8f\n",fun(xk,a,n));
	return 0;
}
