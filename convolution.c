#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "convolution.h"
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

void direct_convolution(PRECISION *x, int nx, PRECISION *h, int nh, PRECISION delta)
{

	PRECISION *x_aux;
	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; //tamano de toda la convolucion
	x_aux = calloc(nx_aux, sizeof(PRECISION));

	int mitad_nh = nh / 2;
   
	//rellenamos el vector auxiliar
	// sobra estamos usando arriba calloc 
	for (k = 0; k < nx_aux; k++)
	{
		x_aux[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		x_aux[k + mitad_nh] = x[k];
	}

	//vamos a tomar solo la convolucion central
	
	for (k = 0; k < nx; k++)
	{
		x[k] = 0;
		
		for (j = 0; j < nh; j++)
		{
			x[k] += h[j] * x_aux[j + k];
		}
		x[k] *= delta;
	}

	free(x_aux);
}

