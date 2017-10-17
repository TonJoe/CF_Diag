#define BEGINNING 500000
#include<iostream>
#include<LA1.h>
#include<time.h>
/**
class CF6
{
public:
	Cdouble J(int n, int m, Cdouble *z);
	Cdouble PJ(qmnumber qm, int nth, Cdouble *z);
	Cdouble CF_Wave(qmnumber *qm, Cdouble *z);
	Cdouble E_Wave(qmnumber *qm, Cdouble *z, Cdouble ze, Cdouble center);

};
class CF5
{
public:
	Cdouble J(int n, int m, Cdouble *z);
	Cdouble PJ(qmnumber qm, int nth, Cdouble *z);
	Cdouble CF_Wave(qmnumber *qm, Cdouble *z);
	Cdouble E_Wave(qmnumber *qm, Cdouble *z, Cdouble ze, Cdouble center);
	Cdouble H_Wave(qmnumber *qm, Cdouble *z, Cdouble center);
};
class CF10
{
public:
	Cdouble J(int n, int m, Cdouble *z);
	Cdouble PJ(qmnumber qm, int nth, Cdouble *z);
	Cdouble CF_Wave(qmnumber *qm, Cdouble *z);
	Cdouble E_Wave(qmnumber *qm, Cdouble *z, Cdouble ze, Cdouble center);
	Cdouble H_Wave(qmnumber *qm, Cdouble *z, Cdouble center);
};**/
class CFN
{
public:
	Cdouble J(int N, int n, int m, Cdouble *z);
	Cdouble PJ(int N, qmnumber qm, int nth, Cdouble *z);
	//Cdouble CFN::LLE(int N, int l, int nth, Cdouble *z); //LLLP Electron with J'
	Cdouble CF_Wave(int N, qmnumber *qm, Cdouble *z);
	Cdouble E_Wave(int N, qmnumber *qm, Cdouble *z, Cdouble ze, Cdouble center);
	Cdouble Landau_Wave1(int N, qmnumber *qm, Cdouble *z, int l);
	Cdouble Landau_Wave2(int N, qmnumber *qm, Cdouble *z, int l);
	//Cdouble H_Wave(int N, qmnumber *qm, Cdouble *z, Cdouble center);
};
