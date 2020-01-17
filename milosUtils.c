#include "lib.h"
#include "defines.h"
#include "milosUtils.h"
//#include "nrutil.h"
//#include "svdcmp.h"
#include "time.h"
#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include <math.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "fftw.h"
#include "convolution.h"
//#include "mkl_vsl.h"
//#include "mkl_lapacke.h"

#define tiempo(ciclos) asm volatile("rdtsc \n\t" \
												: "=A"(ciclos))



extern long long int c1, c2, cd, semi, c1a, c2a, cda; //variables de 64 bits para leer ciclos de reloj
extern long long int c1total, cdtotal;

extern PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
extern int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
extern int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

extern PRECISION *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;

extern PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
extern PRECISION *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
extern PRECISION *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
extern PRECISION *dfi, *dshi;
extern PRECISION *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
extern PRECISION *spectra, *d_spectra, *spectra_mac;
extern PRECISION *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
extern PRECISION *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
extern PRECISION *parcial1, *parcial2, *parcial3;
extern PRECISION *nubB, *nupB, *nurB;
extern PRECISION *G,*GMAC; // VECTOR WITH GAUSSIAN CREATED FOR CONVOLUTION 
extern fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
extern fftw_plan planForwardPSF, planBackwardPSF;

extern fftw_complex * fftw_G_PSF;
//extern VSLConvTaskPtr taskConv;


extern Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original

void spectral_synthesis_convolution(int * nlambda)
{

	int i,j,ishift;
	//int nlambda = NLAMBDA;
	//convolucionamos los perfiles IQUV (spectra)			
	int odd=(*nlambda%2);
	
	int startShift = *nlambda/2;
	if(odd) startShift+=1;

	for (i = 0; i < NPARMS; i++){
		for(j=0;j<* nlambda;j++){
			inSpectraFwPSF[j] = spectra[(*nlambda*i)+j] + 0 * _Complex_I;
		}
		fftw_execute(planForwardPSF);
		// multiplication fft results 
		for(j=0;j<* nlambda;j++){
			inSpectraBwPSF[j] = (outSpectraFwPSF[j]/(*nlambda)) * fftw_G_PSF[j];						
		}
		fftw_execute(planBackwardPSF);
		//shift: -numln/2
		for(j=0,ishift=startShift;j<(*nlambda)/2;j++,ishift++){
			spectra[ishift+i*(*nlambda)]=creal(outSpectraBwPSF[j])*(*nlambda);
		}
		for(j=(*nlambda)/2,ishift=0;j<(*nlambda);j++,ishift++){
			spectra[ishift+i*(*nlambda)]=creal(outSpectraBwPSF[j])*(*nlambda);
		}
	}

}

void response_functions_convolution(int * nlambda)
{

	int i, j, h,k,ishift;
	//int nlambda = NLAMBDA;

	//convolucionamos las funciones respuesta ( d_spectra )

	//PRECISION _Complex *fftaux,*fftaux2, *fftd;
	int odd=(*nlambda%2);
	int startShift = (*nlambda)/2;
	if(odd) startShift+=1;

	for (j = 0; j < NPARMS; j++)
	{
		for (i = 0; i < NTERMS; i++)
		{
			if (i != 7)	{																														 //no convolucionamos S0
				// copy to inSpectra
				for(k=0;k<(*nlambda);k++){
					inSpectraFwPSF[k] = d_spectra[(*nlambda * i + *nlambda * NTERMS * j) + k] + 0 * _Complex_I;
				}
				fftw_execute(planForwardPSF);
				for(h=0;h<(*nlambda);h++){
					inSpectraBwPSF[h] = (outSpectraFwPSF[h]/(*nlambda)) * fftw_G_PSF[h];
				}
				fftw_execute(planBackwardPSF);   			
				//shift 
				for(h=0,ishift=startShift;h<(*nlambda)/2;h++,ishift++){
					d_spectra[ishift+ (*nlambda) * i + (*nlambda) * NTERMS * j]=creal(outSpectraBwPSF[h])*(*nlambda);
				}
				for(h=((*nlambda)/2),ishift=0;h<(*nlambda);h++,ishift++){
					d_spectra[ishift+(*nlambda) * i + (*nlambda) * NTERMS * j]=creal(outSpectraBwPSF[h])*(*nlambda);
				}
			}
		}
	}

}

void AplicaSlight(PRECISION * d_spectra, int numl, PRECISION ALFA, PRECISION * slight){
	int par, il, i;
	// Response Functions 
	for(par=0;par<NPARMS;par++){
		for(il=0;il<NTERMS;il++){
			for(i=0;i<numl;i++){
				d_spectra[numl*il+numl*NTERMS*par+i]=d_spectra[numl*il+numl*NTERMS*par+i]*ALFA;
				if(il==10){ //Magnetic filling factor Response function
					d_spectra[numl*il+numl*NTERMS*par+i]=spectra[numl*par+i]-slight[numl*par+i];
				}
			}
		}

		spectra[i] = spectra[i]*ALFA+slight[i]*(1.0-ALFA);
	}

	for(i=0;i<numl*NPARMS;i++){
		spectra[i] = spectra[i]*ALFA+slight[i]*(1.0-ALFA);
	}
}

void AplicaDelta(Init_Model *model, PRECISION *delta, int *fixed, Init_Model *modelout)
{

	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

	if (fixed[0])
	{
		modelout->eta0 = model->eta0 - delta[0]; // 0
	}
	if (fixed[1])
	{
		if (delta[1] < -300) //300
			delta[1] = -300;
		else if (delta[1] > 300)
			delta[1] = 300;
		modelout->B = model->B - delta[1]; //magnetic field
	}
	if (fixed[2])
	{
		modelout->vlos = model->vlos - delta[2];
	}

	if (fixed[3])
	{
		modelout->dopp = model->dopp - delta[3];
	}

	if (fixed[4])
		modelout->aa = model->aa - delta[4];

	if (fixed[5])
	{
		if (delta[5] < -30) //15
			delta[5] = -30;
		else if (delta[5] > 30)
			delta[5] = 30;

		modelout->gm = model->gm - delta[5]; //5
	}
	if (fixed[6])
	{
		if (delta[6] < -30)
			delta[6] = -30;
		else if (delta[6] > 30)
			delta[6] = 30;

		modelout->az = model->az - delta[6];
	}
	if (fixed[7])
		modelout->S0 = model->S0 - delta[7];
	if (fixed[8])
		modelout->S1 = model->S1 - delta[8];
	if (fixed[9]){
		modelout->mac = model->mac - delta[9]; //9
	}
	if (fixed[10])
		modelout->alfa = model->alfa - delta[10];
}


int check(Init_Model *model)
{
	//Magnetic field
	if (model->B < 0)
	{
		model->B = -(model->B);
		model->gm = 180.0 - (model->gm);
	}
	if (model->B > 4500)
		model->B = 4500;

	//Inclination
	if (model->gm < 0)
		model->gm = -(model->gm);
	if (model->gm > 180)
	{
		model->gm = 360.0 - model->gm;
	}

	//azimuth
	if (model->az < 0)
		model->az = 180 + (model->az); //model->az= 180 + (model->az);
	if (model->az > 180)
	{
		model->az = model->az - 180.0;
	}

	//RANGOS
	//Eta0
	if (model->eta0 < 1)
		model->eta0 = 1;
	if (model->eta0 > 2500) //idl 2500
		model->eta0 = 2500;

	//velocity
	if (model->vlos < (-20)) //20
		model->vlos = (-20);
	if (model->vlos > 20)
		model->vlos = 20;

	//doppler width ;Do NOT CHANGE THIS
	if (model->dopp < 0.0001)
		model->dopp = 0.0001;
	if (model->dopp > 0.6) // idl 0.6
		model->dopp = 0.6;

	if (model->aa < 0.0001) // idl 1e-4
		model->aa = 0.0001;
	if (model->aa > 10.0) //10
		model->aa = 10.0;

	//S0
	if (model->S0 < 0.0001)
		model->S0 = 0.0001;
	if (model->S0 > 10.00)
		model->S0 = 10.00;

	//S1
	if (model->S1 < 0.0001)
		model->S1 = 0.0001;
	if (model->S1 > 10.00)
		model->S1 = 10.00;

	//macroturbulence
	if (model->mac < 0)
		model->mac = 0;
	if (model->mac > 4)
		model->mac = 4;
	
	// filling factor 
	if(model->alfa<0)
		model->alfa = 0.0;
	if(model->alfa>1.0)
		model->alfa = 1.0;
	
	

	return 1;
}


void FijaACeroDerivadasNoNecesarias(PRECISION *d_spectra, int *fixed, int nlambda)
{

	int In, j, i;
	for (In = 0; In < NTERMS; In++)
		if (fixed[In] == 0)
			for (j = 0; j < NPARMS; j++)
				for (i = 0; i < nlambda; i++)
					d_spectra[i + nlambda * In + j * nlambda * NTERMS] = 0;
}


/*
	Tamaño de H es 	 NTERMS x NTERMS
	Tamaño de beta es 1xNTERMS

	return en delta tam 1xNTERMS
*/

int mil_svd(PRECISION *h, PRECISION *beta, PRECISION *delta)
{

	PRECISION epsilon;
	
	static PRECISION h1[NTERMS * NTERMS];
	
	PRECISION *v, *w;
	int i, j;
	//	static PRECISION aux2[NTERMS*NTERMS];
	static PRECISION aux2[NTERMS];
	int aux_nf, aux_nc;
	
	epsilon = 1e-12;
	
	/**/
	
	gsl_vector *eval = gsl_vector_alloc (NTERMS);
  	gsl_matrix *evec = gsl_matrix_alloc (NTERMS, NTERMS);

	for (j = 0; j < NTERMS * NTERMS; j++)
	{
		h1[j] = h[j];
	}

	
	//********************* USING SVD 
	//svdcmp(h1, NTERMS, NTERMS, w, v);

	//******************** USING DIVIDE AND CONQUER
	// TRY USE LAPACKE_dsyevd to calculate eigenvalues and eigenvectors 
	/*MKL_INT n = NTERMS, lda = NTERMS, info;
		// Solve eigenproblem 
	info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR, 'V', 'U', n, h1, lda, w );
		// Check for convergence 
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	*/

	//******************** USING Relatively Robust Representations 
	/*int LDA, LDZ , N;
	LDA = LDZ = N = NTERMS;
	MKL_INT n = N, il, iu, m, lda = LDA, ldz = LDZ, info;
	PRECISION abstol, vl, vu;
	MKL_INT isuppz [2*N];
	PRECISION w3[N], z[LDZ*N];
	abstol = -1.0;
	printf("\n tamaños n %d, m %d , lda %d, ldz %d, abstol %lf \n",n,m,lda,ldz,abstol );
	printf("\n ****\n");
	LAPACKE_dsyevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U', n, h1, lda, vl, vu, il, iu, abstol, &m, w3, z, ldz, isuppz );
	//info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, h1, lda, w );
	printf("\n realizado\n");
	printf("\n**************");*/


	//********************* USING GSL ************************/
	gsl_matrix_view gsl_h1 = gsl_matrix_view_array (h1, NTERMS, NTERMS);
	gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (NTERMS);
	gsl_eigen_symmv(&gsl_h1.matrix, eval, evec, workspace);
	gsl_eigen_symmv_free (workspace);
	w = gsl_vector_ptr(eval,0);
	v = gsl_matrix_ptr(evec,0,0);


	static PRECISION vaux[NTERMS * NTERMS], waux[NTERMS];

	for (j = 0; j < NTERMS * NTERMS; j++)
	{
		vaux[j] = v[j]; 
	}

	for (j = 0; j < NTERMS; j++)
	{
		waux[j] = w[j]; 
	}

	//multmatrixCblas(beta, 1, NTERMS, vaux, NTERMS, NTERMS, aux2, &aux_nf, &aux_nc);
	multmatrix(beta, 1, NTERMS, vaux, NTERMS, NTERMS, aux2, &aux_nf, &aux_nc);

	for (i = 0; i < NTERMS; i++)
	{
      aux2[i]= aux2[i]*((fabs(waux[i]) > epsilon) ? (1/waux[i]): 0.0);
	}

	multmatrix(vaux, NTERMS, NTERMS, aux2, NTERMS, 1, delta, &aux_nf, &aux_nc);
	//multmatrixCblas(vaux, NTERMS, NTERMS, aux2, NTERMS, 1, delta, &aux_nf, &aux_nc);

	
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	return 1;
}



void weights_init(PRECISION *sigma, PRECISION **sigOut, PRECISION noise)
{
	int i;
	//PRECISION *w, *sig;
	PRECISION *sig;

	sig = calloc(4, sizeof(PRECISION));
	if (sigma == NULL)
	{
		for (i = 0; i < 4; i++)
			sig[i] = noise * noise;
	}
	else
	{

		for (i = 0; i < 4; i++)
			sig[i] = (*sigma); // * (*sigma);
	}

	//*wOut = w;
	*sigOut = sig;
}


/*
*
*
* Cálculo de las estimaciones clásicas.
*
*
* lambda_0 :  centro de la línea
* lambda :    vector de muestras
* nlambda :   numero de muesras
* spectro :   vector [I,Q,U,V]
* initModel:  Modelo de atmosfera a ser modificado
*
*
*
* @Author: Juan Pedro Cobos Carrascosa (IAA-CSIC)
*		   jpedro@iaa.es
* @Date:  Nov. 2011
*
*/
void estimacionesClasicas(PRECISION lambda_0, PRECISION *lambda, int nlambda, float *spectro, Init_Model *initModel)
{

	PRECISION x, y, aux, LM_lambda_plus, LM_lambda_minus, Blos, beta_B, Ic, Vlos;
	float *spectroI, *spectroQ, *spectroU, *spectroV;
	PRECISION L, m, gamma, gamma_rad, tan_gamma, C;
	int i;

	//Es necesario crear un lambda en FLOAT para probar como se hace en la FPGA
	PRECISION *lambda_aux = lambda;
	
	

	spectroI = spectro;
	spectroQ = spectro + nlambda;
	spectroU = spectro + nlambda * 2;
	spectroV = spectro + nlambda * 3;

	Ic = spectro[nlambda - 1]; // Continuo ultimo valor de I

	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		aux = (Ic - (spectroI[i] + spectroV[i]));
		x = x + aux * (lambda_aux[i] - lambda_0);
		y = y + aux;
	}

	//Para evitar nan
	if (fabs(y) > 1e-15)
		LM_lambda_plus = x / y;
	else
		LM_lambda_plus = 0;

	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		aux = (Ic - (spectroI[i] - spectroV[i]));
		x = x + aux * (lambda_aux[i] - lambda_0);
		y = y + aux;
	}

	if (fabs(y) > 1e-15)
		LM_lambda_minus = x / y;
	else
		LM_lambda_minus = 0;

	C = (CTE4_6_13 * lambda_0 * lambda_0 * cuantic->GEFF);
	beta_B = 1 / C;

	Blos = beta_B * ((LM_lambda_plus - LM_lambda_minus) / 2);
	Vlos = (VLIGHT / (lambda_0)) * ((LM_lambda_plus + LM_lambda_minus) / 2);


	//------------------------------------------------------------------------------------------------------------
	// //Para probar fórmulación propuesta por D. Orozco (Junio 2017)
	//La formula es la 2.7 que proviene del paper:
	// Diagnostics for spectropolarimetry and magnetography by Jose Carlos del Toro Iniesta and Valent´ýn Mart´ýnez Pillet
	//el 0.08 Es la anchura de la línea en lugar de la resuloción del etalón.

	//Vlos = ( 2*(VLIGHT)*0.08 / (PI*lambda_0)) * atan((spectroI[0]+spectroI[1]-spectroI[3]-spectroI[4])/(spectroI[0]-spectroI[1]-spectroI[3]+spectroI[4]));

	//------------------------------------------------------------------------------------------------------------

	Blos = Blos * 1; //factor de correción x campo debil
	Vlos = Vlos * 1; //factor de correción ...

	//inclinacion
	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		L = fabs(sqrtf(spectroQ[i] * spectroQ[i] + spectroU[i] * spectroU[i]));
		m = fabs((4 * (lambda_aux[i] - lambda_0) * L)); // / (3*C*Blos) ); //2*3*C*Blos mod abril 2016 (en test!)

		x = x + fabs(spectroV[i]) * m;
		y = y + fabs(spectroV[i]) * fabs(spectroV[i]);

	}

	y = y * fabs((3 * C * Blos));

	tan_gamma = fabs(sqrtf(x / y));

	gamma_rad = atan(tan_gamma); //gamma en radianes

	gamma = gamma_rad * (180 / PI); //gamma en grados

	//correccion
	//utilizamos el signo de Blos para ver corregir el cuadrante
	

	if (Blos < 0)
		gamma = (180) - gamma;

	//azimuth

	PRECISION tan2phi, phi;
	int muestra;

	if (nlambda == 6)
		muestra = CLASSICAL_ESTIMATES_SAMPLE_REF;
	else
		muestra = nlambda * 0.75;

	tan2phi = spectroU[muestra] / spectroQ[muestra];

	phi = (atan(tan2phi) * 180 / PI) / 2; //atan con paso a grados

	if (spectroU[muestra] > 0 && spectroQ[muestra] > 0)
		phi = phi;
	else if (spectroU[muestra] < 0 && spectroQ[muestra] > 0)
		phi = phi + 180;
	else if (spectroU[muestra] < 0 && spectroQ[muestra] < 0)
		phi = phi + 90;
	else if (spectroU[muestra] > 0 && spectroQ[muestra] < 0)
		phi = phi + 90;

	PRECISION B_aux;

	B_aux = fabs(Blos / cos(gamma_rad)) * 2; // 2 factor de corrección

	//Vlos = Vlos * 1.5;
	if (Vlos < (-20))
		Vlos = -20;
	if (Vlos > (20))
		Vlos = (20);


	initModel->B = (B_aux > 4000 ? 4000 : B_aux);
	initModel->vlos = Vlos; //(Vlos*1.5);//1.5;
	initModel->gm = gamma;
	initModel->az = phi;
	initModel->S0 = Blos;

	//Liberar memoria del vector de lambda auxiliar
	//free(lambda_aux);
}




/*
 *
 * nwlineas :   numero de lineas espectrales
 * wlines :		lineas spectrales
 * lambda :		wavelength axis in angstrom
			longitud nlambda
 * spectra : IQUV por filas, longitud ny=nlambda
 */

int lm_mils(Cuantic *cuantic, PRECISION *wlines, PRECISION *lambda, int nlambda, float *spectro, int nspectro,
				Init_Model *initModel, PRECISION *spectra, PRECISION *chisqrf,
				PRECISION * slight, PRECISION toplim, int miter, PRECISION *weight, int *fix,
				PRECISION *sigma, PRECISION ilambda, int * INSTRUMENTAL_CONVOLUTION, int * iter)
{

	
	//int iter;
	int i, *fixed, nfree;
	static PRECISION delta[NTERMS];
	
	int iter_chisqr_same;
	PRECISION flambda;
	static PRECISION beta[NTERMS], alpha[NTERMS * NTERMS];
	PRECISION chisqr, ochisqr;
	int clanda, ind;
	Init_Model model;
	//static PRECISION sigmaTemp[4] = {1.0, 1.0, 1.0, 1.0};
	
	nfree = CalculaNfree(nspectro);

	if (nfree == 0)
	{
		return -1; //'NOT ENOUGH POINTS'
	}

	flambda = ilambda;

	if (fix == NULL)
	{
		fixed = calloc(NTERMS, sizeof(PRECISION));
		for (i = 0; i < NTERMS; i++)
		{
			fixed[i] = 1;
		}
	}
	else
	{
		fixed = fix;
	}

	clanda = 0;
	*iter = 0;

	static PRECISION covar[NTERMS * NTERMS];
	static PRECISION betad[NTERMS];

	PRECISION chisqr_mem;

	mil_sinrf(cuantic, initModel, wlines, lambda, nlambda, spectra, AH,slight,spectra_mac, *INSTRUMENTAL_CONVOLUTION);
	me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac, spectra,AH, slight, 0,*INSTRUMENTAL_CONVOLUTION);

	/*mil_sinrf(cuantic, initModel, wlines, lambda, nlambda, spectra, AH,NULL,spectra_mac, 0);
	me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac, AH, slight,0,0);	
	if(*INSTRUMENTAL_CONVOLUTION){
		spectral_synthesis_convolution(&nlambda);
		response_functions_convolution(&nlambda);
	}*/


	/*if(slight!=NULL){
		AplicaSlight(d_spectra,nlambda,initModel->alfa,slight);
	}*/
	//convolucionamos los perfiles IQUV (spectra) y convolucionamos las funciones respuesta ( d_spectra )
	

	/*printf("\nVALORES DE LAS FUNCIONES RESPUESTA MACROTURBULENCIA \n");
	int number_parametros = 0;
	for (number_parametros = 0; number_parametros < NTERMS; number_parametros++)
	{
		if(number_parametros==9){
			printf("\n FUNCION RESPUESTA: %d \n",number_parametros);
			for (int kk = 0; kk < nlambda; kk++)
			{
				printf("1\t%lf\t%le\t%le\t%le\t%le\n", lambda[kk],
				d_spectra[kk + nlambda * number_parametros],
				d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS],
				d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 2],
				d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 3]);
			}
		}
	}
	printf("\n");	

	printf("\n valores de spectro sintetizado INICIAL\n");
	for (int kk = 0; kk < nlambda; kk++)
	{
		printf("1\t%f\t%le\t%le\t%le\t%le\n", lambda[kk], spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
	}*/

	FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);
	covarm(weight, sigma, spectro, nlambda, spectra, d_spectra, beta, alpha);
	//covarm(weight, sigmaTemp, spectro, nlambda, spectra, d_spectra, beta, alpha);

	for (i = 0; i < NTERMS; i++)
		betad[i] = beta[i];

	//printf("\nALPHA: \n");
	for (i = 0; i < NTERMS * NTERMS; i++){
		covar[i] = alpha[i];
		//printf("%f, ",alpha[i]);
	}
	//printf("\n");
	/**************************************************************************/

	ochisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);
	//ochisqr = fchisqr(spectra, nspectro, spectro, weight, sigmaTemp, nfree);

	//printf("\n OBJETIVED CHISQR: %0.10f\n",ochisqr);
	//printf("\n FLAMBDA INICIAL: %lf\n",flambda);
	chisqr_mem = (PRECISION)ochisqr;
	iter_chisqr_same = 0;

	model = *initModel;
	
	do
	{

		for (i = 0; i < NTERMS; i++)
		{
			ind = i * (NTERMS + 1);
			covar[ind] = alpha[ind] * (1.0 + flambda);
		}

		/*printf("\n COVAR \n");
		for (i = 0; i < NTERMS * NTERMS; i++){
			printf("%f, ",covar[i]);
		}
		printf("\n");*/

		mil_svd(covar, betad, delta);
		/*printf("\n DELTA CALCULADO POR EL SVD: \n");
		for( i =0; i< NTERMS;i++){
			printf("%f,",delta[i]);
		}*/
		/*printf("\n BETAD CALCULADO POR EL SVD: \n");
		for( i =0; i< NTERMS;i++){
			printf("%f,",betad[i]);
		}
		printf("\n");*/

		AplicaDelta(initModel, delta, fixed, &model);

		check(&model);

		
		/**************************************************************************/
		//=[Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
		//printf("\n MODELO INICIAL FUNCIONES RESPUESTA: %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",model.eta0,model.B,model.vlos,model.dopp,model.aa,model.gm,model.az,model.S0,model.S1,model.mac,model.alfa);


		//mil_sinrf(cuantic, &model, wlines, lambda, nlambda, spectra, AH,NULL,spectra_mac,0);
		mil_sinrf(cuantic, &model, wlines, lambda, nlambda, spectra, AH,slight,spectra_mac,*INSTRUMENTAL_CONVOLUTION);
		/*printf("\n valores de spectro sintetizado con MIL_SINRF sin aplicar PSF %d\n", *iter);
		for (int kk = 0; kk < nlambda; kk++)
		{
			printf("1\t%f\t%le\t%le\t%le\t%le\n", lambda[kk], spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
		}*/
			
		
		/*if(*INSTRUMENTAL_CONVOLUTION){			
			spectral_synthesis_convolution(&nlambda);
		}*/
		/*printf("\n valores de spectro sintetizado con MIL_SINRF en la iteracion %d\n", *iter);
		for (int kk = 0; kk < nlambda; kk++)
		{
			printf("1\t%f\t%le\t%le\t%le\t%le\n", lambda[kk], spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
		}
		*/
		chisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);
		//chisqr = fchisqr(spectra, nspectro, spectro, weight, sigmaTemp, nfree);
		
		
		/**************************************************************************/
		if(chisqr == chisqr_mem){
			iter_chisqr_same++;
			if(iter_chisqr_same>=3) // exit after 4 iterations with chisqr equal
				clanda = 1;
		}
		else
		{
			chisqr_mem = chisqr;
		}

		//printf("\n CHISQR EN LA ITERACION %d,: %e",*iter,chisqr);
		
		/**************************************************************************/
		if ((fabs((ochisqr-chisqr)*100/chisqr) < toplim) || (chisqr < 0.00001)) // condition to exit of the loop 
			clanda = 1;		
		if (chisqr - ochisqr < 0.)
		{

			flambda = flambda / 10.0;

			*initModel = model;

			//printf("iteration=%d , chisqr = %e CONVERGE	- ilambda= %e \n",*iter,chisqr,flambda);


			//me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac, AH, slight,0,0);
			me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac,spectra, AH, slight,0,*INSTRUMENTAL_CONVOLUTION);

			//convolucionamos las funciones respuesta ( d_spectra )
			
			/*if(*INSTRUMENTAL_CONVOLUTION){
				response_functions_convolution(&nlambda);
			}*/

			/*printf("\nVALORES DE LAS FUNCIONES RESPUESTA MACROTURBULENCIA \n");
			int number_parametros = 0;
			for (number_parametros = 0; number_parametros < NTERMS; number_parametros++)
			{
				if(number_parametros==9){
					printf("\n FUNCION RESPUESTA: %d \n",number_parametros);
					for (int kk = 0; kk < nlambda; kk++)
					{
						printf("1\t%lf\t%le\t%le\t%le\t%le\n", lambda[kk],
						d_spectra[kk + nlambda * number_parametros],
						d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS],
						d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 2],
						d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 3]);
					}
				}
			}
			printf("\n");

			exit(1);*/

			FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);	
			covarm(weight, sigma, spectro, nlambda, spectra, d_spectra, beta, alpha);
			//covarm(weight, sigmaTemp, spectro, nlambda, spectra, d_spectra, beta, alpha);
			

			for (i = 0; i < NTERMS; i++)
				betad[i] = beta[i];

			for (i = 0; i < NTERMS * NTERMS; i++)
				covar[i] = alpha[i];

			ochisqr = chisqr;
			
		}
		else
		{
			flambda = flambda * 10; //10;

			//printf("iteration=%d , chisqr = %e NOT CONVERGE	- lambda= %e \n",*iter,ochisqr,flambda);
		}

		if ((flambda > 1e+7) || (flambda < 1e-25)) 
			clanda=1 ; // condition to exit of the loop 		

		(*iter)++;
		
	} while (*iter < miter && !clanda);

	*chisqrf = ochisqr;
	
	if (fix == NULL)
		free(fixed);

	return 1;
}



/**
 * Make the interpolation between deltaLambda and PSF where deltaLambda es x and PSF f(x)
 *  Return the array with the interpolation. 
 * */
int interpolationSplinePSF(PRECISION *deltaLambda, PRECISION * PSF, PRECISION * lambdasSamples, PRECISION centralLambda, size_t N_PSF, PRECISION * fInterpolated, size_t NSamples){

	size_t i;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
  	gsl_spline *spline_cubic = gsl_spline_alloc(gsl_interp_cspline, N_PSF);
	//gsl_spline *spline_akima = gsl_spline_alloc(gsl_interp_akima, NSamples);
	//gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, NSamples);

	gsl_spline_init(spline_cubic, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_akima, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_steffen, deltaLambda, PSF, N_PSF);

	for (i = 0; i < NSamples; ++i){
   	PRECISION xi = lambdasSamples[i]-centralLambda;
      fInterpolated[i] = gsl_spline_eval(spline_cubic, xi, acc);
      //PRECISION yi_akima = gsl_spline_eval(spline_akima, xi, acc);
      //PRECISION yi_steffen = gsl_spline_eval(spline_steffen, xi, acc);
    }

  gsl_spline_free(spline_cubic);
  //gsl_spline_free(spline_akima);
  //gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);

	return 1;
}


/**
 * Make the interpolation between deltaLambda and PSF where deltaLambda es x and PSF f(x)
 *  Return the array with the interpolation. 
 * */
int interpolationLinearPSF(PRECISION *deltaLambda, PRECISION * PSF, PRECISION * lambdasSamples, PRECISION centralLambda, size_t N_PSF, PRECISION * fInterpolated, size_t NSamples){

	size_t i;
	gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,N_PSF);
   gsl_interp_init(interpolation, deltaLambda, PSF, N_PSF);
   gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();


	for (i = 0; i < NSamples; ++i){
   	PRECISION xi = lambdasSamples[i]-centralLambda;
      fInterpolated[i] = gsl_interp_eval(interpolation, deltaLambda, PSF, xi, accelerator);
    }

  
  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);

	return 1;
}