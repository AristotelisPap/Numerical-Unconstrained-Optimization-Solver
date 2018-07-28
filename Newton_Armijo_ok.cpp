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

double hess_diag_i_i_odd (double xi,double xiplus1,int a)
{
	return (2+8*a*pow(xi,2)-4*a*(xiplus1-xi));
}

double hess_diag_i_i_even (double xi,double xiplus1,int a)
{
	return (2*a);
}

double hess_i_iplus1 (double xi,int a)
{
	return ((-4)*a*xi);
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
		while ((fun(x,a,n)>funpr+theta*lam*fdd)&&(lam>0.005))
		{
			lam=lam/eta;
			for (i=0;i<n;i++)
			{
				x[i]=x[i]-lam*d[i];	
			}
			funx=fun(x,a,n);
		}
	}
	printf("directional=%f,step=%f\n",fdd,lam);
	return (lam);
}



int main()
{
	/* Newton with Armijo and Steepest descent step if necessary*/
	int a,i,j,n,iter=0;
	printf("Choose the parameters a and n:\n");
	scanf("%d %d",&a,&n);
	double x[n],xk[n],dir[n],gradient[n],fun_eval,hessian[n][n];
	double epsilon=pow(10,-5),norm=0,ssgrad=0,step;
	for (i=0;i<n;i++)
	{
		gradient[i]=0;
		for (j=0;j<n;j++)
		{
			hessian[i][j]=0;
		}
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
	printf("initial f(x)=%f\n",fun_eval);
	for (i=0;i<n-1;i=i+2)
	{
		hessian[i][i]=hess_diag_i_i_odd(x[i],x[i+1],a);
		hessian[i][i+1]=hess_i_iplus1(x[i],a);
	}
	for (i=1;i<n;i=i+2)
	{
		hessian[i][i]=hess_diag_i_i_even(x[i],x[i+1],a);
		hessian[i][i-1]=hess_i_iplus1(x[i-1],a);
	}
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("%f\t",hessian[i][j]);
		}
		printf("\n");
	}
	for (i=0;i<n;i=i+2)
	{
		gradient[i]=subgradient_odd(x[i],x[i+1],a);
	}
	for (i=1;i<n;i=i+2)
	{
		gradient[i]=subgradient_even(x[i-1],x[i],a);
	}
	for (i=0;i<n;i++)
	{
		printf("gradient%d=%f\n",i+1,gradient[i]);
	}
	printf("\n");
	for (i=0;i<n;i++)
	{
		ssgrad=ssgrad+pow(gradient[i],2);
	}
	norm=sqrt(ssgrad);
	int s;
	double dchol[n],c[n][n],lchol[n][n];
	while (norm>=epsilon)
	{
		/* Cholesky factorization in order to check p.d.*/
		for (i=0;i<n;i++)
		{
			dchol[i]=0;
			for (j=0;j<n;j++)
			{
				c[i][j]=0;
				lchol[i][j]=0;
			}
		}
		for (i=0;i<n;i++)
		{
			lchol[i][i]=1;
		}
		double sum1,sum2;
		for (j=0;j<n;j++)
		{
			sum1=0;
			for (s=0;s<=j-1;s++)
			{
				sum1=sum1+dchol[s]*pow(lchol[j][s],2);
			}
			c[j][j]=hessian[j][j]-sum1;
			dchol[j]=c[j][j];
			for (i=j+1;i<n;i++)
			{
				sum2=0;
				for (s=0;s<=j-1;s++)
				{
					sum2=sum2+dchol[s]*lchol[i][s]*lchol[j][s];
				}
				c[i][j]=hessian[i][j]-sum2;
				lchol[i][j]=c[i][j]/dchol[j];
			}
		}
		/* Check for p.d. and Steepest Descent step if necessary*/
		int test=0;
		for (i=0;i<n;i++)
		{
			if (dchol[i]<=0)
			{
				test=test+1;
			}
		}
		if (test>0)
		{
			printf("Steepest Descent step\n");
			for (i=0;i<n;i++)
			{
				dir[i]=-gradient[i]/norm;
			}
		}
		else
		/* Solve the system in order to find the direction */
		{
			double delta[n],epsilon[n];
			for (i=0;i<n;i++)
			{
				dir[i]=0;
				delta[i]=0;
				epsilon[i]=0;
			}
			for (i=0;i<n;i++)
			{
				for (j=0;j<=i-1;j++)
				{
					delta[i]=delta[i]-(lchol[i][j]*delta[j]/lchol[i][i]);	
				}
				delta[i]=delta[i]-gradient[i]/lchol[i][i];
			}
			for (i=0;i<n;i++)
			{
				epsilon[i]=delta[i]/dchol[i];
			}
			for (i=n-1;i>=0;i=i-1)
			{
				for (j=i+1;j<n;j++)
				{
					dir[i]=dir[i]-(lchol[j][i]*dir[j]/lchol[i][i]);
				}
				dir[i]=dir[i]+epsilon[i]/lchol[i][i];
			}
		}
		/* Update points */
		step=armijo(x,dir,gradient,a,n);	
		for (i=0;i<n;i++)
		{
			xk[i]=xk[i]+step*dir[i];
			printf("xk%d=%f\n",i+1,xk[i]);
			x[i]=xk[i];
		}
		/* Update gradient and Hessian */
		for (i=0;i<n-1;i=i+2)
		{
			hessian[i][i]=hess_diag_i_i_odd(x[i],x[i+1],a);
			hessian[i][i+1]=hess_i_iplus1(x[i],a);
		}
		for (i=1;i<n;i=i+2)
		{
			hessian[i][i]=hess_diag_i_i_even(x[i],x[i+1],a);
			hessian[i][i-1]=hess_i_iplus1(x[i-1],a);
		}
		for (i=0;i<n;i=i+2)
		{
			gradient[i]=subgradient_odd(x[i],x[i+1],a);
		}
		for (i=1;i<n;i=i+2)
		{
			gradient[i]=subgradient_even(x[i-1],x[i],a);
		}
		ssgrad=0;
		for (i=0;i<n;i++)
		{
			ssgrad=ssgrad+pow(gradient[i],2);
		}
		fun_eval=fun(x,a,n);
		printf("f(x)=%f\n\n",fun_eval);
		norm=sqrt(ssgrad);
		iter=iter+1;
	}
	printf("\tAfter %d iter, we get:\n",iter);
	for (i=0;i<n;i++)
	{
		printf("\tx%d,opt=%2.5f\n",i+1,xk[i]);
	}
	printf("\tf(xopt)=%3.8f\n",fun(xk,a,n));
	return 0;
}
