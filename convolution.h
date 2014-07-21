#ifndef _CONVOLE_H
#define _CONVOLE_H

//#define USE_HEAP

#ifdef USE_HEAP
double __attribute__((optimize("O0"))) *convolution(double *a1, int n1, double *a2, int n2);
#else
double *convolution(double *a1, int n1, double *a2, int n2);
#endif
//double *convolution(double *a1, int n1, double *a2, int n2);
void convolution_free(void);

#endif
