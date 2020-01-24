#include "defines.h"

void direct_convolution(REAL * x, int nx,REAL * h, int nh,PRECISION delta);
void convolve(REAL * Signal, size_t SignalLen, REAL * Kernel, size_t KernelLen);
//int convolutionMKL(double  * h, double * x,  double * y, int  start,VSLConvTaskPtr * task);