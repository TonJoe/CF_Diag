#include<Monca.h>
#include<stdio.h>
Cdouble Metrop(CFN experiment, int N, int steps, qmnumber *qm1, int seed)
{
	Cdouble z[N],r[N];
	Cdouble base1, base1tmp;
	double trial=0.70;
	Cdouble Up(0.0,0.0);
	double niu=1.0/(2.*1.0+1.);		//Need to be ajusted manually.
	double RN=sqrt(2.0*10.0/niu);
	//Cdouble center;
	Cdouble rn(RN,0.0);
	//center=rn*reladistance;
	double count=0;
	srand48(seed);
	for(int i=0;i<N;i++)
	{
		z[i]=polar(drand48()*RN,drand48()*2.*PI);
		r[i]=z[i];
	}
	base1=experiment.CF_Wave(N,qm1,z);
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<N;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<N;i++)
		{
			Cdouble ranmv((drand48()-0.5),(drand48()-0.5));
			//cout<<ranmv<<endl;
			r[i]=r[i]+trial*ranmv;
		}		
		//base1=experiment.CF_Wave(N,qm1,z);	
		base1tmp=experiment.CF_Wave(N,qm1,r);
		if(st>BEGINNING)
		{
		//	Up=Up+conj(base1/base2)*base1/base2;
		//	Up=Up+conj(base1/base2);	
			Up=Up+polar(Vee(N,z),0.0);  //data#1
		//	Up=Up+polar(Vee(6,z),0.0);	//data#2	
		
		}
		if(norm(base1tmp/base1)>drand48())// &&norm(r[st%6])<RN*RN)
		{	
			count=count+1.0;
			base1=base1tmp;
			for(int i=0;i<N;i++)
			{
				z[i]=r[i];
			}
		}
		
	

		if(st%50000==0)
		{
			
			cout<<st<<":"<<count/50000<<"Time: "<<time(NULL)<<endl;
			/**if(count/50000>0.55)
				trial*=1.05;
			else if(count/50000<0.45)
				trial*=0.95;**/
			count=0;
		}
	}
	Up=Up*polar(1./(steps-BEGINNING),0.);	cout<<"Energy="<<Up<<endl;
	
	
	return Up;	

}
Cdouble Metrop(CFN experiment, int N, int steps, qmnumber *qm1, qmnumber *qm2, int l, int seed)
{
	ofstream file;
	file.open("log20_Diag");
	Cdouble z[N],r[N];
	Cdouble base1, base1tmp;
	Cdouble base2, base2tmp;
	double trial=1.0;
	Cdouble Amp1(0.0,0.0),Mo2(0.0,0.0),Overlap21(0.0,0.0);
	double niu=1.0/(2.*1.0+1.);		//Need to be ajusted manually.
	double RN=sqrt(2.0*N/niu);
	Cdouble rn(RN,0.0);
	double count=0;
	srand48(seed);
	for(int i=0;i<N;i++)
	{
		z[i]=polar(drand48()*RN,drand48()*2.*PI);
		r[i]=z[i];
	}
	base1=experiment.CF_Wave(N,qm1,z);
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<N;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<N;i++)
		{
			Cdouble ranmv((drand48()-0.5),(drand48()-0.5));
			r[i]=r[i]+trial*ranmv;
		}
		//Cdouble ranmv((drand48()-0.5),(drand48()-0.5));
		//r[steps%N]=r[steps%N]+trial*ranmv;
		
		base1tmp=experiment.CF_Wave(N,qm1,r);
		base2=experiment.CF_Wave(N,qm2,z);
		if(st>BEGINNING)
		{
		
			Amp1=Amp1+conj(base2/base1)*base2/base1;
			Overlap21=Overlap21+conj(base2/base1);
		}
		if(norm(base1tmp/base1)>drand48())// &&norm(r[st%6])<RN*RN)
		{	
			count=count+1.0;
			base1=base1tmp;
			for(int i=0;i<N;i++)
			{
				z[i]=r[i];
			}
		}
		
	

		if(st%50000==0)
		{
			file<<st<<":"<<count/50000<<"Time: "<<time(NULL)<<endl;
			if(count/50000>0.55)
				trial*=1.005;
			else if(count/50000<0.45)
				trial*=0.995;
			count=0;
		}
	}
	Amp1=Amp1*polar(1./(steps-BEGINNING),0.);
	Overlap21=Overlap21*polar(1./(steps-BEGINNING),0.);
	
	Mo2=Amp1;
	Amp1=sqrt(Amp1);
	Overlap21=Overlap21/Amp1;
	cout<<"Module^2 of |qm2>: "<<Mo2<<" <2|1>: "<<Overlap21<<endl;
	
	file.close();
	return Overlap21;	
}
void Metrop(CFN experiment, int N, int steps, int seed)
{
	qmnumber qm1[N], qm2[N], qm3[N];
	qm1[0].n=2; qm1[0].m=-2;
	qm1[1].n=1; qm1[1].m=-1;
	qm1[2].n=1; qm1[2].m=0;

	qm2[0].n=3; qm2[0].m=-2;
	qm2[1].n=1; qm2[1].m=-1;
	qm2[2].n=1; qm2[2].m=0;

	qm3[0].n=2; qm3[0].m=-2;
	qm3[1].n=2; qm3[1].m=-1;
	qm3[2].n=1; qm3[2].m=0;
	for(int i=3;i<N;i++)
	{
		qm1[i].n=0; qm1[i].m=i-3;
		qm2[i].n=0; qm2[i].m=i-3;
		qm3[i].n=0; qm3[i].m=i-3;
	}

	ofstream file;
	file.open("log_Diag12_2on2on1");
	Cdouble z[N],r[N];
	Cdouble base1, base2, base3, base3tmp;
	double trial=1.0;
	Cdouble Mo13(0.0,0.0),Amp13(0.0,0.0),Amp23(0.0,0.0),Overlap13(0.0,0.0),Overlap23(0.0,0.0);
	double niu=1.0/(2.*1.0+1.);		//Need to be ajusted manually.
	double RN=sqrt(2.0*N/niu);
	Cdouble rn(RN,0.0);
	double count=0;
	srand48(seed);
	for(int i=0;i<N;i++)
	{
		z[i]=polar(drand48()*RN,drand48()*2.*PI);
		r[i]=z[i];
	}
	base3=experiment.CF_Wave(N,qm3,z);

	base1=experiment.CF_Wave(N,qm1,z);
	base2=experiment.CF_Wave(N,qm2,z);
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<N;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<N;i++)
		{
			Cdouble ranmv((drand48()-0.5),(drand48()-0.5));
			r[i]=r[i]+trial*ranmv;
		}		
		base3tmp=experiment.CF_Wave(N,qm3,r);
		
		if(st>BEGINNING)
		{		
			Amp13=Amp13+conj(base1/base3)*base1/base3;
			Overlap13=Overlap13+conj(base1/base3);

			Amp23=Amp23+conj(base2/base3)*base2/base3;
			Overlap23=Overlap23+conj(base2/base3);
		}
		if(norm(base3tmp/base3)>drand48())// &&norm(r[st%6])<RN*RN)
		{	
			count=count+1.0;
			base3=base3tmp;
			for(int i=0;i<N;i++)
			{
				z[i]=r[i];
			}
			base1=experiment.CF_Wave(N,qm1,z);
			base2=experiment.CF_Wave(N,qm2,z);
		}
		
	

		if(st%50000==0)
		{
			file<<st<<":"<<count/50000<<"Time: "<<time(NULL)<<endl;
			if(count/50000>0.55)
				trial*=1.005;
			else if(count/50000<0.45)
				trial*=0.995;
			count=0;
		}
	}
	Amp13=Amp13*polar(1./(steps-BEGINNING),0.);
	Amp23=Amp23*polar(1./(steps-BEGINNING),0.);
	Overlap13=Overlap13*polar(1./(steps-BEGINNING),0.);
	Overlap23=Overlap23*polar(1./(steps-BEGINNING),0.);
	
	Mo13=Amp13;
	Amp13=sqrt(Amp13);
	Amp23=sqrt(Amp23);
	Overlap13=Overlap13/Amp13;
	Overlap23=Overlap23/Amp23;
	cout<<"|<3|3>|^2: "<<Mo13<<" <1|3>: "<<Overlap13<<" <2|3>: "<<Overlap23<<endl;
	
	file.close();	
}
void Module(CFN experiment, int N, int NS, qmnumber **qm, int steps, int seed, Cdouble *Overlap)
{
	ofstream file;
	char filename[10];
	sprintf(filename,"log_%d",N);
	file.open(filename);

	Cdouble *z=new Cdouble[N],*r=new Cdouble[N];
	qmnumber *tempqm=new qmnumber [N];
	double trial=0.5;
	double niu=1.0/(2.*1.0+1.);		//Need to be ajusted manually.
	double RN=sqrt(2.0*N/niu);
	Cdouble rn(RN,0.0);
	double count=0;
	srand48(seed);

	Cdouble *Amp=new Cdouble [NS];
	//Cdouble *Overlap=new Cdouble [NS];
	Cdouble *Module=new Cdouble [NS];
	Cdouble *base=new Cdouble [NS];
	Cdouble base0tmp(0.0,0.0);
	 
	for(int i=0; i<NS; i++)
	{	
		Amp[i]=polar(0.0,0.0);
		Overlap[i]=polar(0.0,0.0);
		base[i]=polar(0.0,0.0);
		Module[i]=polar(0.0,0.0);
	}
	for(int n=0;n<NS;n++)
	{	
		file<<"State #"<<n<<endl;
		for(int i=0;i<N;i++)
		{
			cout<<"**"<<qm[n][i].n<<", "<<qm[n][i].m<<endl;
		}	cout<<endl;
	}

/****Initializing wavefunctions****/
	for(int i=0;i<N;i++)
	{
		z[i]=polar(drand48()*RN,drand48()*2.*PI);
		r[i]=z[i];
	}
	for(int j=0;j<NS;j++)
	{
		//CF_Wave(N,qm1,z);
		base[j]=experiment.CF_Wave(N, qm[j], z);
	}
/****Initialization Done****/
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<N;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<N;i++)
		{
			Cdouble ranmv((drand48()-0.5),(drand48()-0.5));
			r[i]=r[i]+trial*ranmv;
		}	
		base0tmp=experiment.CF_Wave(N,qm[0],r);
		
		if(st>BEGINNING)
		{	
			for(int j=0;j<NS;j++)
			{
				Amp[j]=Amp[j]+conj(base[j]/base[0])*(base[j]/base[0]);
				Overlap[j]=Overlap[j]+conj(base[j]/base[0]);
			}
		}
		if(norm(base0tmp/base[0])>drand48())// &&norm(r[st%6])<RN*RN)
		{	
			count=count+1.0;
			base[0]=base0tmp;
			for(int i=0;i<N;i++)
			{
				z[i]=r[i];
			}
			for(int j=1;j<NS;j++)
			{
				base[j]=experiment.CF_Wave(N, qm[j], z);
			}//cout<<endl;
		}
		
		if(st%50000==0)
		{	
			file<<st<<":"<<count/50000<<"Time: "<<time(NULL)<<endl;
			if(count/50000>0.55)
				trial*=1.005;
			else if(count/50000<0.45)
				trial*=0.995;
			count=0;
		}
	}
	file.close();	//End of log file
	for(int j=0;j<NS;j++)
	{
		Amp[j]=Amp[j]*polar(1./(steps-BEGINNING),0.);
		Overlap[j]=Overlap[j]*polar(1./(steps-BEGINNING),0.);

		Module[j]=Amp[j];
		Amp[j]=sqrt(Amp[j]);
		Overlap[j]=conj(Overlap[j]/Amp[j]);
	}
	cout<<endl;
	
	delete [] tempqm;

	delete [] z;
	delete [] r;
	delete [] Amp;
	//delete [] Overlap;
	delete [] Module;
	delete [] base;
		
}

void Orth(int NS, Cdouble **O, Cdouble **U, Cdouble **NU, Cdouble *Norm)	//Overlap=O_ab=<a|b>, Norm_i=<i'|i'>, U_ab the transform matrix
{
	for(int a=0;a<NS;a++)
	{
		for(int b=0;b<NS;b++)
		{
			if(b==a) U[a][b]=polar(1.0, 0.0);
			else if(b>a) U[a][b]=polar(0.0,0.0);
			else
			{
				U[a][b]=polar(0.0,0.0);
				for(int r=0;r<a;r++)
				{	
					
					Cdouble UO(0.0,0.0);
					for(int d=0;d<=r;d++) UO+=conj(U[r][d])*O[d][a];
					Cdouble UUO(0.0,0.0);
					for(int d=0;d<=r;d++)
						for(int e=0;e<=r;e++) UUO+=conj(U[r][d])*U[r][e]*O[d][e];
					U[a][b]-=UO/UUO*U[r][b];
				}
			}
		}
	}
	for(int a=0;a<NS;a++)
	{
		Norm[a]=polar(0.0,0.0);
		for(int r=0;r<NS;r++)
			for(int d=0;d<NS;d++)
				Norm[a]=Norm[a]+conj(U[a][r])*U[a][d]*O[r][d];
				
	}
	for(int i=0;i<NS;i++)
	{
		for(int j=0;j<NS;j++)
		{
			NU[i][j]=U[i][j]/sqrt(Norm[i]);
		}
	}
}
