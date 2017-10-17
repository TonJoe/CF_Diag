#define NR_END 1
#define LOG2 0.693147180559945
#define PI 3.14159265354
#include <complex>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#define TINY    1.5e-300
#define TINYL    1.5e-1000L

#define FREE_ARG char*
using namespace std;
typedef complex<double> Cdouble;
typedef struct
{
	int n;
	int m;
}qmnumber;
void nrerror(char error_text[]);
void cpxdbl_ludcmp0(complex<double> *a,int n,double *d);
complex<double> cpxdbl_det0(complex<double> *a,int n);
Cdouble Rmv();
double EllipticE(double x);	//Elliptic function, for the use of the evaluation for some wave functions.


/***********************************/
/**Three potentials: e-e, e-b, b-b**/
/***********************************/
double Vee(int n, complex<double> *z);
double Vbb(int n, double niu);
double Vbe(int n, double niu, double RN, complex<double> *z);
Cdouble SingleElectron(Cdouble ze, Cdouble center, double size);
Cdouble Single_Landau(Cdouble ze, int l);
