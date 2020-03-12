#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "convolution.h"
#include "defines.h"
//#include "mkl_vsl.h"

/*

	@autor: Juan Pedro Cobos
	@Date: 31 Marzo 2011
	@loc: IAA- CSIC
	
	Convolucion para el caso Sophi: convolucion central de x con h.
	
	direct_convolution(x,h,delta)
	x spectro
	h gaussiana -perfil instrumental- � ojo, solo con longitud impar!
	delta anchura de muestreo

	--nota: como h es simetrico no se invierte su orden 
	
	//result=calloc(nresult,sizeof(PRECISION*));
	//free();

	_juanp
*/
extern PRECISION *dirConvPar;
//extern REAL *resultConv;

void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh)
{

	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central

	for (k = 0; k < nx; k++)
	{
		//x[k] = 0;
		double aux = 0;
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
		}
		x[k] = aux;
	}
}

void direct_convolution(REAL *x, int nx, PRECISION *h, int nh)
{

	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central
	double aux;
	for (k = 0; k < nx; k++)
	{
		//x[k] = 0;
		aux = 0;
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
		}
		x[k] = aux;
	}
}

void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh, REAL *result, int delta)
{

	int nx_aux;
	int k, j;
	
	nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central

	for (k = 0; k < nx; k++)
	{

		double aux = 0;
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
		}
		result[k] = aux * delta;
	}
}

void convolve(REAL *Signal, size_t SignalLen, double *Kernel, size_t KernelLen, REAL *Result, int delta)
{
	size_t n;
	size_t kmin, kmax, k;
	double aux;
	size_t SignalLen_Aux = SignalLen + SignalLen;
	size_t KernelLen_Aux = KernelLen + KernelLen;
	REAL *signalAux = calloc(SignalLen_Aux, sizeof(REAL));
	double *kernelAux = calloc(KernelLen_Aux, sizeof(double));

	size_t dirConvParLen = (KernelLen_Aux + KernelLen_Aux - 1);
	double * dirConvParAux = calloc(dirConvParLen,sizeof(double));
	int mitad_nh = SignalLen / 2;
	for (n = mitad_nh; n < SignalLen + mitad_nh; n++)
	{
		signalAux[n] = Signal[n-mitad_nh];
		kernelAux[n] = Kernel[n-mitad_nh];
	}
	printf("\n result convolution: ");
	for (n = 0; n < SignalLen_Aux + KernelLen_Aux - 1; n++){
		//dirConvPar[n] = 0;
		aux = 0;
		kmin = (n >= KernelLen_Aux - 1) ? n - (KernelLen_Aux - 1) : 0;
		kmax = (n < SignalLen_Aux - 1) ? n : SignalLen_Aux - 1;

		for (k = kmin; k <= kmax; k++){
			aux += signalAux[k] * kernelAux[n - k];
		}
		dirConvParAux[n] = aux;
		printf("\n %e",aux);
	}
	printf("\n");
	/*REAL *signalAux = calloc(SignalLen + SignalLen - 1, sizeof(REAL));
	double *ketnelAux = calloc(KernelLen + KernelLen - 1, sizeof(double));
	int mitad_nh = SignalLen / 2;
	for (n = mitad_nh; n < SignalLen + mitad_nh; n++)
	{
		signalAux[n] = Signal[n];
		ketnelAux[n] = Kernel[n];
	}
	for (n = 0; n < SignalLen + KernelLen - 1; n++)
	{
		dirConvPar[n] = 0;
		for (i = 0; i <= n; i++)
		{
			dirConvPar[n] = dirConvPar[n] + (signalAux[i] * ketnelAux[n - i]);
		}
	}
	*/
	
	for (n = 0; n < SignalLen; n++){
		Result[n] = delta*dirConvParAux[n + mitad_nh+SignalLen];
	}
}

/**
 * Method to do circular convolution over signal 'x'. We assume signal 'x' and 'h' has the same size. 
 * The result is stored in array 'result'
 * 
 * 
 * */

void convCircular(REAL *x, double *h, int size, REAL *result)
{
	int i,j,ishift,mod;
	double aux;

	int odd=(size%2);		
	int startShift = size/2;
	if(odd) startShift+=1;	
	ishift = startShift;

	for(i=0; i < size ; i++){
		aux = 0;
    	for(j=0; j < size; j++){
			mod = i-j;
			if(mod<0)
				mod = size+mod;
			aux += h[j] * x[mod];
    	}
		//resultConv[i] = aux;
		if(i < size/2)
			result[ishift++] = aux;
		else
			result[i-(size/2)] = aux;		
	}
	/*for(i=0;i<size;i++){
		result[i] = resultConv[i];
	}*/

}



/*void convCircular(REAL *x, double *h, int size, REAL *result)
{
	int i,j, k,ishift,mod;
	double aux;

	int odd=(size%2);		
	int startShift = size/2;
	if(odd) startShift+=1;	
	ishift = startShift;

	for(i=0; i < size ; i++){
		aux = 0;
    	for(j=0; j < size; j++){
			mod = i-j;
			if(mod<0)
				mod = size+mod;
			aux += h[j] * x[mod];
    	}
		//resultConv[i] = aux;
		if(i < size/2)
			resultConv[ishift++] = aux;
		else
			resultConv[i-(size/2)] = aux;		
	}
	for(i=0;i<size;i++){
		result[i] = resultConv[i];
	}*/

	/*for(i=0,ishift=startShift;i<n/2;i++,ishift++){
		result[ishift]=resultConv[i];
	}
	for(i=(n/2),ishift=0;i<n;i++,ishift++){
		result[ishift]=resultConv[i];
	} */ 

	
   /*dirConvPar[0]=h[0];

   for(j=1;j<n;j++)            //folding h(n) to h(-n)
   	dirConvPar[j]=h[n-j];

    //Circular convolution
	aux = 0;
	for(i=0;i<n;i++)
   	aux+=x[i]*dirConvPar[i];

	resultConv[ishift++] = aux;
	//resultConv[0] = aux;

	for(k=1;k<n;k++){
   	
      //circular shift

      for(j=1;j<n;j++)
      	dirConvPar2[j]=dirConvPar[j-1];
		
      dirConvPar2[0]=dirConvPar[n-1];
		aux=0;
      for(i=0;i<n;i++){
      	dirConvPar[i]=dirConvPar2[i];
         aux+=x[i]*dirConvPar2[i];
      }
		//resultConv[k] = aux;
		if(k <m/2)
			resultConv[ishift++] = aux;
		else
			resultConv[k-(n/2)] = aux;
	}
	for(i=0;i<n;i++){
		result[i] = resultConv[i];
	}	
	*/
	

//}
