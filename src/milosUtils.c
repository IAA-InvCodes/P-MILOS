#include "lib.h"
#include "defines.h"
#include "milosUtils.h"
#include "time.h"
#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include <math.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include "convolution.h"



extern PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
extern int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
extern int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

extern REAL *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;

extern REAL *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
extern REAL *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
extern REAL *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
extern REAL *dfi, *dshi;
extern REAL *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
//extern PRECISION *spectra, *d_spectra, *spectra_mac;

extern REAL *spectra, *d_spectra, *spectra_mac;
extern REAL *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
extern REAL *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
extern REAL *parcial1, *parcial2, *parcial3;
extern REAL *nubB, *nupB, *nurB;
extern PRECISION *G,*GMAC;
//extern REAL *G;
extern fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
extern fftw_plan planForwardPSF, planBackwardPSF;

extern fftw_complex * fftw_G_PSF;
//extern VSLConvTaskPtr taskConv;


extern Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original


extern gsl_vector *eval;
extern gsl_matrix *evec;
extern gsl_eigen_symmv_workspace * workspace;

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

void AplicaSlight(REAL * d_spectra, int numl, PRECISION ALFA, PRECISION * slight){
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

	if (fixed[0])  // ETHA 0 
	{
		modelout->eta0 = model->eta0 - delta[0]; // 0
	}
	if (fixed[1]) // B
	{
		if (delta[1] < -300) //300
			delta[1] = -300;
		else if (delta[1] > 300)
			delta[1] = 300;
		modelout->B = model->B - delta[1]; //magnetic field
	}
	if (fixed[2]) // VLOS
	{
		/*if (delta[2] > 2)
			delta[2] = 2;

		if (delta[2] < -2)
			delta[2] = -2;		*/
		modelout->vlos = model->vlos - delta[2];
	}

	if (fixed[3]) // DOPPLER WIDTH
	{

		/*if (delta[3] > 1e-2)
			delta[3] = 1e-2;
		else if (delta[3] < -1e-2)
			delta[3] = -1e-2;*/
		modelout->dopp = model->dopp - delta[3];
	}

	if (fixed[4]) // DAMPING 
		modelout->aa = model->aa - delta[4];

	if (fixed[5])  // GAMMA 
	{
		if (delta[5] < -30) //15
			delta[5] = -30;
		else if (delta[5] > 30)
			delta[5] = 30;

		modelout->gm = model->gm - delta[5]; //5
	}
	if (fixed[6]) // AZIMUTH
	{
		/*if (delta[6] < -15)
			delta[6] = -15;
		else if (delta[6] > 15)
			delta[6] = 15;*/

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

	// damping 
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


void FijaACeroDerivadasNoNecesarias(REAL *d_spectra, int *fixed, int nlambda)
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
	static PRECISION aux2[NTERMS];
	int aux_nf, aux_nc;
	
	epsilon = 1e-12;

	for (j = 0; j < NTERMS * NTERMS; j++)
	{
		h1[j] = h[j];
	}

	gsl_matrix_view gsl_h1 = gsl_matrix_view_array (h1, NTERMS, NTERMS);
	gsl_eigen_symmv(&gsl_h1.matrix, eval, evec, workspace);
	w = gsl_vector_ptr(eval,0);
	v = gsl_matrix_ptr(evec,0,0);

	multmatrix(beta, 1, NTERMS, v, NTERMS, NTERMS, aux2, &aux_nf, &aux_nc);

	for (i = 0; i < NTERMS; i++)
	{
		aux2[i]= aux2[i]*((fabs(w[i]) > epsilon) ? (1/w[i]): 0.0);
	}

	multmatrix(v, NTERMS, NTERMS, aux2, NTERMS, 1, delta, &aux_nf, &aux_nc);
	
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
void estimacionesClasicas(PRECISION lambda_0, PRECISION *lambda, int nlambda, float *spectro, Init_Model *initModel, int forInitialUse)
{

	double x, y, aux, LM_lambda_plus, LM_lambda_minus, Blos, beta_B, Ic, Icmax, Vlos;
	double aux_vlos,x_vlos,y_vlos;
	float *spectroI, *spectroQ, *spectroU, *spectroV;
	double L, m, gamma, gamma_rad, tan_gamma, C;
	int i;

	spectroI = spectro;
	spectroQ = spectro + nlambda;
	spectroU = spectro + nlambda * 2;
	spectroV = spectro + nlambda * 3;


	Ic = spectro[nlambda - 1]; // Continuo ultimo valor de I

	/*Icmax = spectro[0];
	int index =0;
	for (i = 0; i < nlambda; i++)
	{
		if(spectroI[i]>Ic){
			Icmax = spectroI[i];
			index = i;
		}
	}*/

	x = 0;
	y = 0;
	x_vlos = 0;
	y_vlos = 0;
	for (i = 0; i < nlambda-1 ; i++)
	{
		aux = (Ic - (spectroI[i] + spectroV[i]));
		aux_vlos = (Ic - spectroI[i]);
		x += (aux * (lambda[i] - lambda_0));
		x_vlos += (aux_vlos * (lambda[i] - lambda_0));
		y += aux;
		y_vlos += aux_vlos;
	}

	//Para evitar nan
	if (fabs(y) > 1e-15)
		LM_lambda_plus = x / y;
	else
		LM_lambda_plus = 0;

	x = 0;
	y = 0;
	for (i = 0; i < nlambda-1 ; i++)
	{
		aux = (Ic - (spectroI[i] - spectroV[i]));
		x += (aux * (lambda[i] - lambda_0));
		y += aux;
	}

	if (fabs(y) > 1e-15)
		LM_lambda_minus = x / y;
	else
		LM_lambda_minus = 0;

	C = (CTE4_6_13 * (lambda_0*lambda_0) * cuantic->GEFF);
	beta_B = 1 / C;

	Blos = (1 / C) * ((LM_lambda_plus - LM_lambda_minus) / 2);
	//Vlos = (VLIGHT / (lambda_0)) * ((LM_lambda_plus + LM_lambda_minus) / 2);
	Vlos = (VLIGHT / (lambda_0)) * ((x_vlos/y_vlos) / 2); // for now use the center without spectroV only spectroI 


	//------------------------------------------------------------------------------------------------------------
	// //Para probar fórmulación propuesta por D. Orozco (Junio 2017)
	//La formula es la 2.7 que proviene del paper:
	// Diagnostics for spectropolarimetry and magnetography by Jose Carlos del Toro Iniesta and Valent´ýn Mart´ýnez Pillet
	//el 0.08 Es la anchura de la línea en lugar de la resuloción del etalón.

	//Vlos = ( 2*(VLIGHT)*0.08 / (PI*lambda_0)) * atan((spectroI[0]+spectroI[1]-spectroI[3]-spectroI[4])/(spectroI[0]-spectroI[1]-spectroI[3]+spectroI[4]));

	//------------------------------------------------------------------------------------------------------------

	//inclinacion
	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		L = FABS(SQRT(spectroQ[i] * spectroQ[i] + spectroU[i] * spectroU[i]));
		m = fabs((4 * (lambda[i] - lambda_0) * L)); // / (3*C*Blos) ); //2*3*C*Blos mod abril 2016 (en test!)

		x = x + FABS(spectroV[i]) * m;
		y = y + FABS(spectroV[i]) * FABS(spectroV[i]);

	}

	y = y * fabs((3 * C * Blos));

	tan_gamma = fabs(sqrt(x / y));

	gamma_rad = atan(tan_gamma); //gamma en radianes
	gamma = gamma_rad * (180 / PI); //gamma en grados

	if(forInitialUse){
		if(gamma>=85 && gamma <=90){  
			gamma_rad = 85 *(PI/180);
		}
		if(gamma>90 && gamma <=95){ 
			gamma_rad = 95 *(PI/180);
		}
	}
	//correction 
	//we use the sign of Blos to see to correct the quadrant
	if (Blos < 0)
		gamma = (180) - gamma;


	// CALCULATIONS FOR AZIMUTH 
	PRECISION tan2phi, phi;

	double sum_u =0.0, sum_q = 0.0;
	for(i=0;i<nlambda;i++){
		if( fabs(spectroU[i]>0.0001 || fabs(spectroQ[i])>0.0001  )){
			sum_u += spectroU[i];
			sum_q += spectroQ[i];
		}
	}
	tan2phi = sum_u/sum_q;
	phi = (atan(tan2phi) * 180 / PI) / 2;
	if ( sum_u > 0 && sum_q > 0 )
		phi = phi;
	else if ( sum_u < 0 && sum_q > 0 )
		phi = phi + 180;
	else if ( sum_u< 0 && sum_q < 0 )
		phi = phi + 90;
	else if ( sum_u > 0 && sum_q < 0 )
		phi = phi + 90;
	
	// END CALCULATIONS FOR AZIMUTH 
	
	PRECISION B_aux;
	B_aux = fabs(Blos / cos(gamma_rad)); // 

	
	if (Vlos < (-4))
		Vlos = -4;
	if (Vlos > (4))
		Vlos = 4;


	initModel->B = (B_aux > 4000 ? 4000 : B_aux);
	initModel->vlos = Vlos;
	initModel->gm = gamma;
	initModel->az = phi;

	if(!forInitialUse) // store Blos in SO if we are in non-initialization use
		initModel->S0 = Blos;

	//Liberar memoria del vector de lambda auxiliar
	
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
				Init_Model *initModel, REAL *spectra, float *chisqrf,
				PRECISION * slight, PRECISION toplim, int miter, REAL *weight, int *fix,
				REAL *sigma, REAL ilambda, int * INSTRUMENTAL_CONVOLUTION, int * iter)
{

	

	REAL PARBETA_better = 5.0;
   	REAL PARBETA_worst = 10.0;
  	REAL PARBETA_FACTOR = 1.0;

	//int iter;
	int i, *fixed, nfree;
	static PRECISION delta[NTERMS];
	
	
	REAL flambda;
	static REAL beta[NTERMS], alpha[NTERMS * NTERMS];
	REAL chisqr, ochisqr, chisqr0;
	int clanda, ind;
	Init_Model model;	
	
	nfree = CalculaNfree(nspectro);

	if (nfree == 0)
	{
		return -1; //'NOT ENOUGH POINTS'
	}

	flambda = ilambda;

	if (fix == NULL)
	{
		static int fixed_static[NTERMS];
		fixed  = fixed_static;
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

	

	mil_sinrf(cuantic, initModel, wlines, lambda, nlambda, spectra, AH,slight,spectra_mac, *INSTRUMENTAL_CONVOLUTION);
	me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac, spectra,AH, slight, *INSTRUMENTAL_CONVOLUTION);

	FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);
	covarm(weight, sigma, spectro, nlambda, spectra, d_spectra, beta, alpha);

	for (i = 0; i < NTERMS; i++)
		betad[i] = beta[i];

	for (i = 0; i < NTERMS * NTERMS; i++){
		covar[i] = alpha[i];
	}


	ochisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);
	chisqr0 = ochisqr;

	model = *initModel;
	
	do
	{
		
		// CHANGE VALUES OF DIAGONAL 
		for (i = 0; i < NTERMS; i++)
		{
			ind = i * (NTERMS + 1);
			covar[ind] = alpha[ind] * (1.0 + flambda);
		}

		mil_svd(covar, betad, delta);
		AplicaDelta(initModel, delta, fixed, &model);

		check(&model);
		mil_sinrf(cuantic, &model, wlines, lambda, nlambda, spectra, AH,slight,spectra_mac,*INSTRUMENTAL_CONVOLUTION);
	
		chisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);
		
		/**************************************************************************/

		//printf("\n CHISQR EN LA ITERACION %d,: %e",*iter,chisqr);
		
		/**************************************************************************/
		if ((FABS((ochisqr-chisqr)*100/chisqr) < toplim) || (chisqr < 0.0001)) // condition to exit of the loop 
			clanda = 1;		
		if (chisqr - ochisqr < 0.)
		{

			//flambda = flambda / 10.0;
			flambda=flambda/(PARBETA_better*PARBETA_FACTOR);
			*initModel = model;
			me_der(cuantic, initModel, wlines, lambda, nlambda, d_spectra, spectra_mac,spectra, AH, slight,*INSTRUMENTAL_CONVOLUTION);
			FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);	
			covarm(weight, sigma, spectro, nlambda, spectra, d_spectra, beta, alpha);
			
			for (i = 0; i < NTERMS; i++)
				betad[i] = beta[i];

			for (i = 0; i < NTERMS * NTERMS; i++)
				covar[i] = alpha[i];

			ochisqr = chisqr;
			
		}
		else
		{
			//flambda = flambda * 10; //10;
			flambda=flambda*PARBETA_worst*PARBETA_FACTOR;
		}

		if ((flambda > 1e+7) || (flambda < 1e-25)) 
			clanda=1 ; // condition to exit of the loop 		

		(*iter)++;
		PARBETA_FACTOR = log10f(chisqr)/log10f(chisqr0);

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
int interpolationSplinePSF(PRECISION *deltaLambda, PRECISION * PSF, PRECISION * lambdasSamples, size_t N_PSF, PRECISION * fInterpolated, size_t NSamples){

	size_t i;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
  	gsl_spline *spline_cubic = gsl_spline_alloc(gsl_interp_cspline, N_PSF);
	//gsl_spline *spline_akima = gsl_spline_alloc(gsl_interp_akima, NSamples);
	//gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, NSamples);

	gsl_spline_init(spline_cubic, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_akima, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_steffen, deltaLambda, PSF, N_PSF);

	for (i = 0; i < NSamples; ++i){
   	
      //fInterpolated[i] = gsl_spline_eval(spline_cubic, xi, acc);
      //PRECISION yi_akima = gsl_spline_eval(spline_akima, xi, acc);
      //PRECISION yi_steffen = gsl_spline_eval(spline_steffen, lambdasSamples[i], acc);
		PRECISION yi = gsl_spline_eval(spline_cubic, lambdasSamples[i], acc);
		if(!gsl_isnan(yi)){
			fInterpolated[i] = yi;
		}
		else
		{
			fInterpolated[i] = 0.0f;
		}
		
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
int interpolationLinearPSF(PRECISION *deltaLambda, PRECISION * PSF, PRECISION * lambdasSamples, size_t N_PSF, PRECISION * fInterpolated, size_t NSamples,double offset){

	size_t i;
	gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,N_PSF);
   	gsl_interp_init(interpolation, deltaLambda, PSF, N_PSF);
   	gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

	for (i = 0; i < NSamples; ++i){
		//printf("\n VALOR A INERPOLAR EN X %f, iteration %i\n",lambdasSamples[i],i);
		double aux;
		if(offset>0){
			if(lambdasSamples[i]-offset>= deltaLambda[0]){
				aux = gsl_interp_eval(interpolation, deltaLambda, PSF, lambdasSamples[i]-offset, accelerator);
						// if lambdasSamples[i] is out of range from deltaLambda then aux is GSL_NAN, we put nan values to 0. 
				if(!gsl_isnan(aux)) 
					fInterpolated[i] = aux;
				else
					fInterpolated[i] = 0.0f;
			}
			else
			{
				aux = 0.0f;
			}
		}
		else{
			if(lambdasSamples[i]+offset<= deltaLambda[NSamples-1]){
				aux = gsl_interp_eval(interpolation, deltaLambda, PSF, lambdasSamples[i]+offset, accelerator);
						// if lambdasSamples[i] is out of range from deltaLambda then aux is GSL_NAN, we put nan values to 0. 
				if(!gsl_isnan(aux)) 
					fInterpolated[i] = aux;
				else
					fInterpolated[i] = 0.0f;
			}
			else
			{
				aux = 0.0f;
			}			
		}
   }

  	// normalizations 
	double cte = 0;
	for(i=0; i< NSamples; i++){
		cte += fInterpolated[i];
	}
	for(i=0; i< NSamples; i++){
		fInterpolated[i] /= cte;
	}
  	gsl_interp_accel_free(accelerator);
	gsl_interp_free(interpolation);
  	
	/*for(i=0; i< NSamples; i++){
		if(fInterpolated[i]<1e-3)
			fInterpolated[i] =0;
	}*/


	return 1;
}
