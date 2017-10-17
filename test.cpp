#define _NUMBER_ELECTRON 14
#define _NUMBER_STATE 14

#include<Monca.h>
#include <fstream>
using namespace std;
int main()
{
	qmnumber qm1[_NUMBER_ELECTRON], qm2[_NUMBER_ELECTRON];
		
	qmnumber **qm=new qmnumber*[_NUMBER_STATE];
	for(int i=0;i<_NUMBER_STATE;i++)
		qm[i]=new qmnumber[_NUMBER_ELECTRON];
	
	qm[0][0].n=2; qm[0][0].m=-2;
	qm[0][1].n=1; qm[0][1].m=-1;
	qm[0][2].n=1; qm[0][2].m=0;
	for(int i=3;i<_NUMBER_ELECTRON;i++)
	{
		qm[0][i].n=0; qm[0][i].m=i-3;
	}

	for(int j=1;j<=3;j++)
	{
		qm[j][0].n=2; qm[j][0].m=-2;
		qm[j][1].n=1; qm[j][1].m=-1;
		qm[j][2].n=1; qm[j][2].m=0;
		for(int i=3;i<_NUMBER_ELECTRON;i++)
		{
			qm[j][i].n=0; qm[j][i].m=i-3;
		}
		qm[j][j-1].n=qm[j][j-1].n+1;

	}
	for(int j=4;j<_NUMBER_STATE;j++)
	{
		qm[j][0].n=2; qm[j][0].m=-2;
		qm[j][1].n=1; qm[j][1].m=-1;
		qm[j][2].n=1; qm[j][2].m=0;
		for(int i=3;i<_NUMBER_ELECTRON;i++)
		{
			qm[j][i].n=0; qm[j][i].m=i-3;
		}
		qm[j][j].n=qm[j][j].n+1;

	}
	CFN experiment;
	//Cdouble Overlaps[_NUMBER_STATE];
	//Cdouble **coef;
	Cdouble **Overlaps=new Cdouble*[_NUMBER_STATE];
	for(int i=0;i<_NUMBER_STATE;i++)
		Overlaps[i]=new Cdouble[_NUMBER_STATE];
	Cdouble **U=new Cdouble*[_NUMBER_STATE];
	for(int i=0;i<_NUMBER_STATE;i++)
		U[i]=new Cdouble[_NUMBER_STATE];
	Cdouble *Norm=new Cdouble[_NUMBER_STATE];
	
	Cdouble **NU=new Cdouble*[_NUMBER_STATE];
	for(int i=0;i<_NUMBER_STATE;i++)
		NU[i]=new Cdouble[_NUMBER_STATE];
	
for(int times=0; times<3;times++)
{
	for(int j=0;j<_NUMBER_STATE;j++)
		Module(experiment, _NUMBER_ELECTRON, _NUMBER_STATE-j, qm+j, 5000000, times, Overlaps[j]+j);	//Overlap matrix only have upper half. need to make up for the lower half.
	for(int i=0;i<_NUMBER_STATE;i++)
		for(int j=0;j<i;j++)
			Overlaps[i][j]=conj(Overlaps[j][i]);
	
	for(int i=0;i<_NUMBER_STATE;i++)
	{	
		cout<<endl;
		for(int j=0;j<_NUMBER_STATE;j++)
			cout<<Overlaps[i][j]<<" ";
	}cout<<endl;
	
	Orth(_NUMBER_STATE, Overlaps, U, NU, Norm);
	for(int i=0;i<_NUMBER_STATE;i++)
	{	
		cout<<endl;
		for(int j=0;j<_NUMBER_STATE;j++)
			cout<<U[i][j]<<" ";
	}cout<<endl<<endl;
	for(int i=0;i<_NUMBER_STATE;i++)
	{
		cout<<Norm[i]<<"  ";
	}cout<<endl;

	for(int i=0;i<_NUMBER_STATE;i++)
	{	
		cout<<"{";
		for(int j=0;j<_NUMBER_STATE;j++)
		{
			cout<<"("<<NU[i][j].real()<<"+I*"<<NU[i][j].imag()<<") ";
			if(j<_NUMBER_STATE-1)
			cout<<",";
		}
		cout<<"}"<<endl;
	}
}
	for(int i=0;i<_NUMBER_STATE;i++)
		delete [] qm[i];
	delete [] qm;
	for(int i=0;i<_NUMBER_STATE;i++)
		delete [] Overlaps[i];
	delete [] Overlaps;
	for(int i=0;i<_NUMBER_STATE;i++)
		delete [] U[i];
	delete [] U;
	for(int i=0;i<_NUMBER_STATE;i++)
		delete [] NU[i];
	delete [] NU;
 
	delete [] Norm;
}
