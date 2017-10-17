#include<LA1.h>
using namespace std;
void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run_time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
void cpxdbl_ludcmp0(complex<double> *a,int n,double *d)
{
	double veevee[n+1];
	int i,imax,j,k;
	double big,dum,temp;
	complex<double> sum,dum2;

  *d=1.0;
  for(i=0;i<n;i++)
  {
    big=0.0;
    for(j=0;j<n;j++)
      if( (temp=abs(a[i*n+j]))>big)       big=temp;
    if(big==0.0)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0/big;      /* Save the scaling */
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum= a[i*n+j];
      for(k=0;k<i;k++)     sum=sum-a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++){
      sum=a[i*n+j];
      for(k=0;k<j;k++)
        sum=sum-a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*abs(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=0;k<n;k++){
        dum2=a[imax*n+k];
        a[imax*n+k]=a[j*n+k];
  a[imax*n+k]=a[j*n+k];
        a[j*n+k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }
    if(abs(a[j*n+j])==0.0)     a[j*n+j]=TINY+polar(0.0,1.0)*TINY;
    if(j!=n){
      dum2=1.0/a[j*n+j];
      for(i=j+1;i<n;i++)      a[i*n+j]=a[i*n+j]*dum2;
    }
  }
}

complex<double> cpxdbl_det0(complex<double> *a,int n)
{
  int i;
  double d;
  complex<double> det;

  cpxdbl_ludcmp0(a,n,&d);
  det=d;
  for(i=0;i<n;i++){
    det=det*a[i*n+i];
  }
return det;
}
Cdouble Rmv()
{
	double r=drand48();
	if(r<0.25)
	{
		Cdouble c(1.,1.);
		return c;
	}
	else if(0.25<=r&&r<0.5)
	{
		Cdouble c(1.,-1.);
		return c;
	}	
	else if(0.5<=r&&r<0.75)
	{
		Cdouble c(-1.,1.);
		return c;
	}
	else
	{
		Cdouble c(-1.,-1.);
		return c;
	}
}
/***********************************/
/**Three potentials: e-e, e-b, b-b**/
/***********************************/
double Vee(int n, complex<double> *z)
{
	double V=0;
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++)
		{
			V=V+1./abs(z[i]-z[j]);
		}
	return V;
}
double Vbb(int n, double niu)
{
	return n*8./3./PI*sqrt(niu*n/2);
}
double Vbe(int n, double niu, double RN, complex<double> *z)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum=sum+EllipticE((abs(z[i])/RN)*(abs(z[i])/RN));
	}
	return -sqrt(2.*niu*n)*sum*2/PI;
}

double EllipticE(double x)	//Elliptic function, for the use of the evaluation for some wave functions.
{
	return PI/2.-PI/8.*x-3.*PI/128.*x*x-5.*PI/512.*pow(x,3)-175.*PI/32768.*pow(x,4)-441.*PI/131072.*pow(x,5)-4851*PI/2097152.*pow(x,6)-14157.*PI/8388608.*pow(x,7)-2760615.*PI/2147483648.*pow(x,8)-8690825.*PI/8589934592*pow(x,9)-112285459.*PI/137438953472.*pow(x,10);

}
Cdouble SingleElectron(Cdouble ze, Cdouble center, double size)
{	
	Cdouble result(exp(-0.25*norm(ze-center)),0.0);
	
	double Rx=real(center), Ry=imag(center);
	double zx=real(ze),zy=-imag(ze);
	Cdouble exponent(0.0,0.5*(Rx*zy-Ry*zx));
	exponent=exp(exponent);
	return result*exponent; 
	//return exponent;
}
Cdouble Single_Landau(Cdouble ze, int l)
{
	return pow(ze, l)*exp(-0.25*norm(ze));
}
