#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include <math.h>
#include <stdio.h>
#include "fftw.h"
#include "defines.h"
/*
	
	direc : FFT_FORWARD (-1) or FFT_BACKWARD (+1)
*/

PRECISION _Complex *fft_d(PRECISION *spectra, int nspectra, int direc)
{

	fftw_complex *in, *out;
	int i;

	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nspectra);

	for (i = 0; i < nspectra; i++)
	{
		in[i] = spectra[i] + 0 * _Complex_I;
	}

	out = fft_c(in, nspectra, direc);

	fftw_free(in);

	return (PRECISION _Complex *) out;
}

PRECISION _Complex *fft_c(PRECISION _Complex *spectra, int nspectra, int direc)
{

	fftw_complex *out;
	fftw_plan p;
	int i;

	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nspectra);

	p = fftw_plan_dft_1d(nspectra, spectra, out, direc, FFTW_ESTIMATE);

	fftw_execute(p);

	if (direc == FFT_FORWARD)
	{
		for (i = 0; i < nspectra; i++)
		{
			out[i] = out[i] / nspectra;
		}
	}

	fftw_destroy_plan(p);

	return (PRECISION _Complex *)out;
}



//Convout's array size must be 2*size.
void FFTConvolve( PRECISION * data,  PRECISION * kernel, int n, int h)
{
    
    int i;

    //Pad the 2nd half of the data with 0's.
    int sizeTot = n + h -1;
	 int mitad_nh = h / 2;
	 fftw_complex * in_sequence, * freq_sequence, * in_data, * freq_data, * rev_data, * time_data;

    in_sequence = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);
	 freq_sequence = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);
	 in_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);
	 freq_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);
	 rev_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);
	 time_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * sizeTot);

    fftw_plan p1 = fftw_plan_dft_1d(sizeTot, in_sequence, freq_sequence, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan p2 = fftw_plan_dft_1d(sizeTot, in_data, freq_data,         FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan rev = fftw_plan_dft_1d(sizeTot, rev_data, time_data,       FFTW_BACKWARD, FFTW_ESTIMATE);

    for( i = 0; i < n; i++ )
    {
        
        in_data[i] = data[i] + 0 * _Complex_I;
    }
    for( ; i < sizeTot; i++ )
    {
        in_data[i] = 0 ;
        
    }

    for( i = 0; i < h; i++ )
    {
        in_sequence[i] = kernel[i] + 0 * _Complex_I;
        
    }
    for( ; i < sizeTot; i++ )
    {
        in_sequence[i] = 0 ;
        
    }

    fftw_execute(p1);
    fftw_execute(p2);

    for( i = 0; i < sizeTot; i++ )
    {  
		  rev_data[i] = (freq_data[i]*freq_sequence[i])/sizeTot;
    }

    fftw_execute(rev);

    for( i = 0; i < n; i++ )
    {
        data[i] = creal(time_data[i+mitad_nh]); // The i term should be non-existent.
    }

	 fftw_free(in_sequence);
	 fftw_free(freq_sequence);
	 fftw_free(in_data);
	 fftw_free(freq_data);
	 fftw_free(rev_data);
	 fftw_free(time_data);
	 fftw_destroy_plan(p1);
	 fftw_destroy_plan(p2);
	 fftw_destroy_plan(rev);
}




void FFTConvolve1D( PRECISION * data, int sizeData, int h, fftw_complex * kernel, fftw_complex * in_data, fftw_complex * freq_data, fftw_complex * rev_data, fftw_complex * time_data, fftw_plan p2, fftw_plan rev)
{
    
    int i;

    //Pad the 2nd half of the data with 0's.
    int sizeTot = sizeData + h -1;
	 int mitad_nh = h / 2;

    for( i = 0; i < sizeData; i++ )
    {        
        in_data[i] = data[i] + 0 * _Complex_I;
    }
    for( ; i < sizeTot; i++ )
    {   
        in_data[i] = 0 ;
    }

    
    fftw_execute(p2);

    for( i = 0; i < sizeTot; i++ )
    {
		  rev_data[i] = (freq_data[i]*kernel[i])/sizeTot;
    }

    fftw_execute(rev);

    for( i = 0; i < sizeData; i++ )
    {
		 //printf("\n%lf",creal(time_data[i]));
        data[i] = creal(time_data[i+mitad_nh]); // The i term should be non-existent.
    }
}