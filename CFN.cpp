#include<CF.h>
Cdouble CFN::J(int N, int n, int m, Cdouble *z)	//N: particle number. n: lambda level.
{
	//Cdouble neg1(1.0, PI);
	//double Norm;
	//Norm=1.0/sqrt(4.0*pow(2.0,double(m))*pow(double(n+1.0),double(m)));
	//cout<<Norm;
	Cdouble product(1.0,0.0),product1(0.0,0.0);//,product2(0.0,0.0);
	for(int i=1;i<N;i++)
	{
		product=product*(z[0]-z[i]);
	}
	if(n==0)
	{	
			
		return pow(z[0],m+n)*product;
	}
	if(n==1)
	{	
		for(int i=1;i<N;i++)
			product1=product1+1.0/(z[0]-z[i]);
		return -pow(z[0],m+n)*product1*product;
	}
	if(n==2)
	{	

		for(int i=1;i<N;i++)
			for(int j=i+1;j<N;j++)
				product1=product1+2.0/((z[0]-z[i])*(z[0]-z[j]));
		return pow(z[0],m+n)*product1*product;
	}
	if(n==3)
	{	

		for(int i=1;i<N;i++)
			for(int j=i+1;j<N;j++)
				for(int k=j+1;k<N;k++)
					product1=product1+6.0/((z[0]-z[i])*(z[0]-z[j])*(z[0]-z[k]));
		return -pow(z[0],m+n)*product1*product;
	}
}
/**New CF_Wave function
Cdouble CFN::J(int N, int n, int m, Cdouble *z)	//N: particle number. n: lambda level.
{
	Cdouble product1(0.0,0.0);//,product2(0.0,0.0);
	if(n==0)
	{				
		return pow(z[0],m+n);
	}
	if(n==1)
	{	
		for(int i=1;i<N;i++)
			product1=product1+1.0/(z[0]-z[i]);
		return -pow(z[0],m+n)*product1;
	}
	if(n==2)
	{	

		for(int i=1;i<N;i++)
			for(int j=i+1;j<N;j++)
				product1=product1+2.0/((z[0]-z[i])*(z[0]-z[j]));
		return pow(z[0],m+n)*product1;
	}
	if(n==3)
	{	

		for(int i=1;i<N;i++)
			for(int j=i+1;j<N;j++)
				for(int k=j+1;k<N;k++)
					product1=product1+6.0/((z[0]-z[i])*(z[0]-z[j])*(z[0]-z[k]));
		return -pow(z[0],m+n)*product1;
	}
}**/
Cdouble CFN::PJ(int N, qmnumber qm, int nth, Cdouble *z)
{
	Cdouble r[N];
	for(int i=0;i<N;i++)
	{
		r[i]=z[(i+nth)%N];
		
	}//cout<<"***************"<<qm.n<<"  "<<qm.m<<endl;
	return J(N, qm.n, qm.m, r);
}

Cdouble CFN::CF_Wave(int N, qmnumber *qm, Cdouble *z)
{
	double Zoom=exp(-0.5*N*sqrt(N));
	Cdouble matrix[N*N];
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			matrix[i*N+j]=PJ(N, qm[i],j,z)*Zoom;
		}
	}
	double sumnorm=0.0;
	for(int i=0;i<N;i++)
	{
		sumnorm=sumnorm+norm(z[i]);
	}
	sumnorm=exp(-0.25*sumnorm);
	Cdouble expo(sumnorm,0.0);
	Cdouble result;
	result=cpxdbl_det0(matrix, N)*expo;
	return result;
}
/**
Cdouble CFN::E_Wave(int N, qmnumber *qm, Cdouble *z, Cdouble ze, Cdouble center)	//
{
	Cdouble r[N+1],sum(0.0,0.0);
	sum=SingleElectron(ze,center,1.0)*CF_Wave(qm, z); //cout<<SingleElectron(ze,center,1.0);
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			r[j]=z[j];
		}
		r[N]=r[i];
		r[i]=ze;
		sum=sum-SingleElectron(r[N],center,1.0)*CF_Wave(qm, r);
	}
	return sum;
}
**/

Cdouble CFN::Landau_Wave1(int N, qmnumber *qm, Cdouble *z, int l)	//N-1 CFs & 1 e-
{
	
	Cdouble r[N],sum(0.0,0.0),tmp;
	sum=Single_Landau(z[N-1],l)*CF_Wave(N-1,qm, z); //cout<<SingleElectron(ze,center,1.0);
	for(int i=0;i<N-1;i++)
	{
		for(int j=0;j<N;j++)
		{
			r[j]=z[j];
		}
		tmp=r[N-1];
		r[N-1]=r[i];
		r[i]=tmp;
		sum=sum-Single_Landau(r[N-1],l)*CF_Wave(N-1, qm, r);
	}
	return sum;
}
Cdouble CFN::Landau_Wave2(int N, qmnumber *qm, Cdouble *z, int l)	//N-1 CFs & 1 e-
{
	

	double Zoom=exp(-0.5*N*sqrt(N));
	Cdouble product(1.0,0.0);
	Cdouble matrix[N*N];
	for(int j=0;j<N;j++)
	{
		product=polar(1.0,0.0);
		for(int k=1;k<N;k++)
			product=product*(z[(j+k)%N]-z[j]);
		matrix[j]=Single_Landau(z[j],l)/product;
	}
	for(int i=1;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			matrix[i*N+j]=PJ(N, qm[i-1],j,z)*Zoom;
		}
	}
	double sumnorm=0.0;
	for(int i=0;i<N;i++)
	{
		sumnorm=sumnorm+norm(z[i]);
	}
	sumnorm=exp(-0.25*sumnorm);
	Cdouble expo(sumnorm,0.0);
	Cdouble result;
	result=cpxdbl_det0(matrix, N)*expo;
	return result;
}
