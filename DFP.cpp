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
	double value=0;
	value=(-4)*a*xi*(xiplus1-pow(xi,2))+2*(xi-1);
    return (value);
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
	/* DFP with Armijo using equation (6.15) p.139 Nocedal & Wright*/
	int a,i,j,n,iter=0;
	printf("Choose the parameters a and n:\n");
	scanf("%d %d",&a,&n);
	double x[n],xk[n],dir[n],gradient[n],gradientk[n],fun_eval,H_dfp[n][n];
	double epsilon=pow(10,-5),norm=0,ssgrad=0,step,s[n],y[n],test;
	double Hy[n],ytH[n],HyytH[n][n],ytHy,frac1[n][n],sst[n][n],frac2[n][n],ssdir;
	for (i=0;i<n;i++)
	{
		s[i]=0;
		y[i]=0;
		dir[i]=0;
		gradientk[i]=0;
		gradient[i]=0;
		for (j=0;j<n;j++)
		{
			H_dfp[i][j]=0;
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
	/* Take the initial H_dfp to be identity */
	for (i=0;i<n;i++)
	{
		H_dfp[i][i]=1;
	}
	printf("\nInitialization\n");
	printf("f(x)=%f\n",fun(x,a,n));
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("%f\t",H_dfp[i][j]);
		}
		printf("\n");
	}
	fun_eval=fun(x,a,n);
	for (i=0;i<n;i=i+2)
	{
		gradient[i]=subgradient_odd(x[i],x[i+1],a);
		gradientk[i]=gradient[i];
	}
	for (i=1;i<n;i=i+2)
	{
		gradient[i]=subgradient_even(x[i-1],x[i],a);
		gradientk[i]=gradient[i];
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
	while (norm>=epsilon)
	{
		/* Find initial direction */
		for (i=0;i<n;i++)
		{
			dir[i]=0;
			for (j=0;j<n;j++)
			{
				dir[i]=dir[i]-H_dfp[i][j]*gradient[j];
			}
		}
		ssdir=0;
		for (i=0;i<n;i++)
		{
			ssdir=ssdir+pow(dir[i],2);
		}
		ssdir=sqrt(ssdir);
		for (i=0;i<n;i++)
		{
			dir[i]=dir[i]/ssdir;
		}
		for (i=0;i<n;i++)
		{
			printf("dir%d=%f\n",i+1,dir[i]);
		}
		step=armijo(x,dir,gradient,a,n);
		/* Initial update of points */
		for (i=0;i<n;i++)
		{
			x[i]=xk[i];
			x[i]=x[i]+step*dir[i];
		}
		/* Define sk and yk and check if their inner product is positive */
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
			s[i]=x[i]-xk[i];
			y[i]=gradient[i]-gradientk[i];		
		}
		test=0;
		for (i=0;i<n;i++)
		{
			test=test+s[i]*y[i];
		}
		/* DFP Inverse Hessian Update */
		if (test>0)
		{
			/* Create Hk*yk */
			for (i=0;i<n;i++)
			{
				Hy[i]=0;
				for (j=0;j<n;j++)
				{
					Hy[i]=Hy[i]+H_dfp[i][j]*y[j];
				}
			}
			/* Create yk(transpose)*Hk */
			for (i=0;i<n;i++)
			{
				ytH[i]=0;
				for (j=0;j<n;j++)
				{
					ytH[i]=ytH[i]+y[j]*H_dfp[j][i];
				}
			}
			/* Create Hk*yk*yk(transpose)*Hk */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					HyytH[i][j]=Hy[i]*ytH[j];
				}
			}
			/* Create yk(transpose)*Hk*yk */
			ytHy=0;
			for (i=0;i<n;i++)
			{
				ytHy=ytHy+ytH[i]*y[i];
			}
			/* Create the 1st fraction of the equation 6.15 */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					frac1[i][j]=HyytH[i][j]/ytHy;
				}
			}
			/* Create sk*sk(transpose) */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					sst[i][j]=s[i]*s[j];
				}
			}
			/* Create the 2nd fraction of the equation 6.15 */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					frac2[i][j]=sst[i][j]/test;
				}
			}
			/* Update Hk+1 */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					H_dfp[i][j]=H_dfp[i][j]-frac1[i][j]+frac2[i][j];
				}
			}
		}
		else
		{
			/* if sk(transpose)*yk<=0 then do not do the update! */
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					H_dfp[i][j]=H_dfp[i][j];
				}
			}
		}
		for (i=0;i<n;i++)
		{
			xk[i]=x[i];
			printf("xk%d=%f\n",i+1,xk[i]);
			gradientk[i]=gradient[i];
		}
		ssgrad=0;
		for (i=0;i<n;i++)
		{
			ssgrad=ssgrad+pow(gradient[i],2);
		}
		fun_eval=fun(x,a,n);
		printf("f(x)=%f\n",fun_eval);
		norm=sqrt(ssgrad);
		iter=iter+1;
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
