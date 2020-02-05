#include "defines.h"

void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh);
void direct_convolution(REAL *x, int nx, PRECISION *h, int nh);
void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh,REAL * result,int delta);
void convolve(REAL * Signal, size_t SignalLen, double * Kernel, size_t KernelLen, REAL * Result , int delta);
void convCircular(REAL *x, double *h, int size, REAL *result);
