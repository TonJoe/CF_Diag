#include<CF.h>
#include <fstream>
Cdouble Metrop(CFN experiment, int N, int steps, qmnumber *qm1, int seed);
Cdouble Metrop(CFN experiment, int N, int steps, qmnumber *qm1, qmnumber *qm2, int l, int seed);
void Metrop(CFN experiment, int N, int steps, int seed);
void Module(CFN experiment, int N, int NS, qmnumber **qm, int steps, int seed, Cdouble *Overlaps);
void Orth(int P_N, Cdouble **Overlaps, Cdouble **U, Cdouble **NU, Cdouble *Norm);
