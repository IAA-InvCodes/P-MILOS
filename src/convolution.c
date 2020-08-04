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
extern REAL *resultConv;

void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh)
{

	//int nx_aux;
	int k, j;

	//nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	/*for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}*/

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

	//int nx_aux;
	int k, j;

	//nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	/*for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}*/

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


void direct_convolution_ic(REAL *x, int nx, PRECISION *h, int nh, REAL Ic)
{

	//int nx_aux;
	int k, j;

	//nx_aux = nx + nh - 1; // tamano de toda la convolucion

	int mitad_nh = nh / 2;

	// rellenamos el vector auxiliar
	/*for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}*/

	/*for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = Ic - x[k];
	}*/

	// vamos a tomar solo la convolucion central
	double aux;
	for (k = 0; k < nx; k++)
	{
		//x[k] = 0;
		aux = 0;
		int N_start_point=k-(nh/2);
		for (j = 0; j < nh; j++)
		{
			if(((N_start_point+j)>=0) && ((N_start_point+j)<nh))
				aux += h[j] * (Ic-x[N_start_point+j]);
			//aux += h[j] * dirConvPar[j + k];
		}
		x[k] = Ic - aux;
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
		if(i < size/2)
			resultConv[ishift++] = aux;
		else
			resultConv[i-(size/2)] = aux;
	}
	for(i=0;i<size;i++){
		result[i] = resultConv[i];
	}

}

