
//    _______             _______ _________ _        _______  _______
//   (  ____ \           (       )\__   __/( \      (  ___  )(  ____ \
//   | (    \/           | () () |   ) (   | (      | (   ) || (    \/
//   | |         _____   | || || |   | |   | |      | |   | || (_____
//   | |        (_____)  | |(_)| |   | |   | |      | |   | |(_____  )
//   | |                 | |   | |   | |   | |      | |   | |      ) |
//   | (____/\           | )   ( |___) (___| (____/\| (___) |/\____) |
//   (_______/           |/     \|\_______/(_______/(_______)\_______)
//
//
// CMILOS v1.0 (2020)
// RTE INVERSION C code (based on the ILD code MILOS by D. Orozco)
// Manuel (IAA-CSIC)
//

/*
;      eta0 = line-to-continuum absorption coefficient ratio
;      B = magnetic field strength       [Gauss]
;      vlos = line-of-sight velocity     [km/s]
;      dopp = Doppler width              [Angstroms]
;      aa = damping parameter
;      gm = magnetic field inclination   [deg]
;      az = magnetic field azimuth       [deg]
;      S0 = source function constant
;      S1 = source function gradient
;      mac = macroturbulent velocity     [km/s]
;      alpha = filling factor of the magnetic component [0->1]

*/

#include <time.h>
#include "defines.h"
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "milosUtils.h"
#include "lib.h"
#include "readConfig.h"
#include <unistd.h>
#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original


PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

REAL *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
REAL *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
//PRECISION gp4_gp2_rhoq[NLAMBDA],gp5_gp2_rhou[NLAMBDA],gp6_gp2_rhov[NLAMBDA];
REAL *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;
REAL *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
REAL *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
REAL *dfi, *dshi;
REAL CC, CC_2, sin_gm, azi_2, sinis, cosis, cosis_2, cosi, sina, cosa, sinda, cosda, sindi, cosdi, sinis_cosa, sinis_sina;
REAL *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
REAL *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
REAL *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
REAL *parcial1, *parcial2, *parcial3;
REAL *nubB, *nupB, *nurB;
REAL **uuGlobalInicial;
REAL **HGlobalInicial;
REAL **FGlobalInicial;


PRECISION *GMAC,* GMAC_DERIV;
PRECISION * dirConvPar;
//REAL *resultConv;
PRECISION * G = NULL;
//REAL * G;

REAL AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];

REAL * opa;
int FGlobal, HGlobal, uuGlobal;

//PRECISION *d_spectra, *spectra, *spectra_mac;
REAL *d_spectra, *spectra, *spectra_mac;



// GLOBAL variables to use for FFT calculation 
//VSLConvTaskPtr taskConv;
fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
fftw_complex * inSpectraFwMAC, *inSpectraBwMAC, *outSpectraFwMAC, *outSpectraBwMAC;
fftw_plan planForwardPSF, planBackwardPSF;
fftw_plan planForwardMAC, planBackwardMAC;
fftw_complex * inFilterMAC, * inFilterMAC_DERIV, * outFilterMAC, * outFilterMAC_DERIV;
fftw_plan planFilterMAC, planFilterMAC_DERIV;


fftw_complex * fftw_G_PSF, * fftw_G_MAC_PSF, * fftw_G_MAC_DERIV_PSF;
fftw_complex * inPSF_MAC, * inMulMacPSF, * inPSF_MAC_DERIV, *inMulMacPSFDeriv, *outConvFilters, * outConvFiltersDeriv;
fftw_plan planForwardPSF_MAC, planForwardPSF_MAC_DERIV,planBackwardPSF_MAC, planBackwardPSF_MAC_DERIV;


//Convolutions values
int sizeG = 0;
PRECISION FWHM = 0;

ConfigControl configCrontrolFile;

// fvoigt memory consuption
 _Complex double *z,* zden, * zdiv;


// 
gsl_vector *eval;
gsl_matrix *evec;
gsl_eigen_symmv_workspace * workspace;

int main(int argc, char **argv)
{
	int i; // for indexes
	PRECISION *wlines;
	int nlambda;
	Init_Model *vModels;
	float chisqrf, * vChisqrf;
	int * vNumIter; // to store the number of iterations used to converge for each pixel
	int indexLine; // index to identify central line to read it 

	
	// allocate memory for eigen values
	eval = gsl_vector_alloc (NTERMS);
  	evec = gsl_matrix_alloc (NTERMS, NTERMS);
	workspace = gsl_eigen_symmv_alloc (NTERMS);

	//*****
	Init_Model INITIAL_MODEL;
	PRECISION * deltaLambda, * PSF;
	PRECISION initialLambda, step, finalLambda;
	int N_SAMPLES_PSF;
	int posWL=0;
	//----------------------------------------------

	PRECISION * slight = NULL;
	int dimStrayLight;

	const char  * nameInputFileSpectra ;
	
	
	char nameOutputFilePerfiles [4096];
	
	const char	* nameInputFileLines;
	
	const char	* nameInputFilePSF ;	

   FitsImage * fitsImage;
	PRECISION  dat[7];

	/********************* Read data input from file ******************************/

	/* Read data input from file */

	loadInitialValues(&configCrontrolFile);
	readTrolFile(argv[1],&configCrontrolFile,1);

	nameInputFileSpectra = configCrontrolFile.ObservedProfiles;
	nameInputFileLines = configCrontrolFile.AtomicParametersFile;
	
	nameInputFilePSF = configCrontrolFile.PSFFile;
	FWHM = configCrontrolFile.FWHM;


	/***************** READ INIT MODEL ********************************/
	if(configCrontrolFile.InitialGuessModel[0]!='\0' && !readInitialModel(&INITIAL_MODEL,configCrontrolFile.InitialGuessModel)){
		printf("\nERROR READING GUESS MODEL 1 FILE\n");
		exit(EXIT_FAILURE);
	}
	
	/***************** READ WAVELENGHT FROM GRID OR FITS ********************************/
	PRECISION * vLambda, *vOffsetsLambda;

	if(configCrontrolFile.useMallaGrid){ // read lambda from grid file
      indexLine = readMallaGrid(configCrontrolFile.MallaGrid, &initialLambda, &step, &finalLambda, 1);      
      nlambda = ((finalLambda-initialLambda)/step)+1;
		vOffsetsLambda = calloc(nlambda,sizeof(PRECISION));
		vOffsetsLambda[0] = initialLambda;
		for(i=1;i<nlambda;i++){
			vOffsetsLambda[i] = vOffsetsLambda[i-1]+step;
		}
      // pass to armstrong 
      initialLambda = initialLambda/1000.0;
      step = step/1000.0;
      finalLambda = finalLambda/1000.0;
	   vLambda = calloc(nlambda,sizeof(PRECISION));
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(nameInputFileLines,dat,indexLine,1);
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n QUANTUM LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
			exit(1);
		}
		vLambda[0]=configCrontrolFile.CentralWaveLenght+(initialLambda);
   	for(i=1;i<nlambda;i++){
        	vLambda[i]=vLambda[i-1]+step;
     	}
	}
	else{ // read lambda from fits file
		vLambda = readFitsLambdaToArray(configCrontrolFile.WavelengthFile,&indexLine,&nlambda);
		if(vLambda==NULL){
			printf("\n FILE WITH WAVELENGHT HAS NOT BEEN READ PROPERLY, please check it.\n");
			free(vLambda);
			exit(EXIT_FAILURE);
		}
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(nameInputFileLines,dat,indexLine,1);
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n QUANTUM LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
			exit(1);
		}		
	}

	/*********************************************** INITIALIZE VARIABLES  *********************************/

	CC = PI / 180.0;
	CC_2 = CC * 2;

	wlines = (PRECISION *)calloc(2, sizeof(PRECISION));
	wlines[0] = 1;
	wlines[1] = configCrontrolFile.CentralWaveLenght;

	/******************* CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/

	cuantic = create_cuantic(dat,1);

	/****************************************************************************************************/
	int numln=nlambda;
	// MACROTURBULENCE PLANS
	inFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	outFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	planFilterMAC = fftw_plan_dft_1d(numln, inFilterMAC, outFilterMAC, FFT_FORWARD, FFTW_MEASURE );
	inFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	outFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	planFilterMAC_DERIV = fftw_plan_dft_1d(numln, inFilterMAC_DERIV, outFilterMAC_DERIV, FFT_FORWARD, FFTW_MEASURE );


	inSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	outSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	planForwardMAC = fftw_plan_dft_1d(numln, inSpectraFwMAC, outSpectraFwMAC, FFT_FORWARD, FFTW_MEASURE );
	inSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
	outSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);		
	planBackwardMAC = fftw_plan_dft_1d(numln, inSpectraBwMAC, outSpectraBwMAC, FFT_BACKWARD, FFTW_MEASURE );

	// ********************************************* IF PSF HAS BEEN SELECTEC IN TROL READ PSF FILE OR CREATE GAUSSIAN FILTER ***********//
	if(configCrontrolFile.ConvolveWithPSF){
		
		if(configCrontrolFile.FWHM > 0){
			//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
			G = fgauss_WL(FWHM,vLambda[1]-vLambda[0],vLambda[0],vLambda[nlambda/2],nlambda,&sizeG);
		}else{
			// read the number of lines 
			FILE *fp;
			char ch;
			N_SAMPLES_PSF=0;
			//open file in read more
			fp=fopen(nameInputFilePSF,"r");
			if(fp==NULL)
			{
				printf("File \"%s\" does not exist!!!\n",nameInputFilePSF);
				//slog_error(0,"File \"%s\" does not exist!!!\n",nameInputFilePSF);
				return 0;
			}

			//read character by character and check for new line	
			while((ch=fgetc(fp))!=EOF)
			{
				if(ch=='\n')
					N_SAMPLES_PSF++;
			}
			
			//close the file
			fclose(fp);
			if(N_SAMPLES_PSF>0){
				deltaLambda = calloc(N_SAMPLES_PSF,sizeof(PRECISION));
				PSF = calloc(N_SAMPLES_PSF,sizeof(PRECISION));
				readPSFFile(deltaLambda,PSF,nameInputFilePSF,configCrontrolFile.CentralWaveLenght);
				// CHECK if values of deltaLambda are in the same range of vLambda. For do that we truncate to 4 decimal places 
				if( (trunc(vOffsetsLambda[0])) < (trunc(deltaLambda[0]))  || (trunc(vOffsetsLambda[nlambda-1])) > (trunc(deltaLambda[N_SAMPLES_PSF-1])) ){
					printf("\n\n ERROR: The wavelength range given in the PSF file is smaller than the range in the mesh file [%lf,%lf] [%lf,%lf]  \n\n",deltaLambda[0],vOffsetsLambda[0],deltaLambda[N_SAMPLES_PSF-1],vOffsetsLambda[nlambda-1]);
					exit(EXIT_FAILURE);
				}
				G = malloc(nlambda * sizeof(PRECISION));
				
				double offset=0;
				for(i=0;i<nlambda && !posWL;i++){
					//if( (trunc(vLambda[i]*1000)/1000)== (trunc(configCrontrolFile.CentralWaveLenght*1000)/1000))
					if( fabs(trunc(vOffsetsLambda[i]))==0) 
						posWL = i;
				}
				if(posWL!= (nlambda/2)){ // move center to the middle of samples
					//printf("\nPOS CENTRAL WL %i",posWL);
					offset = (((nlambda/2)-posWL)*step)*1000;
					//printf ("\n OFFSET IS %f\n",offset);
				}
				
				interpolationLinearPSF(deltaLambda,  PSF, vOffsetsLambda ,N_SAMPLES_PSF, G, nlambda,offset);

				sizeG = nlambda;
			}
			else{
				//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
				//PRECISION * fgauss_WL(PRECISION FWHM, PRECISION step_between_lw, PRECISION lambda0, PRECISION lambdaCentral, int nLambda, int * sizeG)
				//G = fgauss_WL(0.75,vLambda[1]-vLambda[0],vLambda[0],vLambda[nlambda/2],nlambda,&sizeG);
				printf("\n****************** ERROR THE PSF FILE is empty or damaged.******************\n");
				exit(EXIT_FAILURE);
			}

		}
		
		//PSF FILTER PLANS 
		inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planForwardPSF = fftw_plan_dft_1d(numln, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_MEASURE );
		inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);		
		planBackwardPSF = fftw_plan_dft_1d(numln, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_MEASURE );

		fftw_complex * in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numln);
		int i;
		for (i = 0; i < numln; i++)
		{
			in[i] = G[i] + 0 * _Complex_I;
		}
		fftw_G_PSF = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numln);
		fftw_plan p = fftw_plan_dft_1d(numln, in, fftw_G_PSF, FFT_FORWARD, FFTW_MEASURE );
		fftw_execute(p);
		for (i = 0; i < numln; i++)
		{
			fftw_G_PSF[i] = fftw_G_PSF[i] / numln;
		}
		fftw_destroy_plan(p);
		fftw_free(in);
		
		inPSF_MAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		fftw_G_MAC_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planForwardPSF_MAC = fftw_plan_dft_1d(numln, inPSF_MAC, fftw_G_MAC_PSF, FFT_FORWARD, FFTW_MEASURE );
		inMulMacPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outConvFilters = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planBackwardPSF_MAC = fftw_plan_dft_1d(numln, inMulMacPSF, outConvFilters, FFT_BACKWARD, FFTW_MEASURE );


		inPSF_MAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		fftw_G_MAC_DERIV_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planForwardPSF_MAC_DERIV = fftw_plan_dft_1d(numln, inPSF_MAC_DERIV, fftw_G_MAC_DERIV_PSF, FFT_FORWARD, FFTW_MEASURE );
		inMulMacPSFDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outConvFiltersDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planBackwardPSF_MAC_DERIV = fftw_plan_dft_1d(numln, inMulMacPSFDeriv, outConvFiltersDeriv, FFT_BACKWARD, FFTW_MEASURE );			

	}		



	/****************************************************************************************************/
	//  IF NUMBER OF CYCLES IS LES THAN 0 THEN --> WE USE CLASSICAL ESTIMATES 
	//  IF NUMBER OF CYCLES IS 0 THEN -->  DO SYNTHESIS FROM THE INIT MODEL 
	//  IF NUMBER OF CYCLES IS GREATER THAN 0 --> READ FITS FILE OR PER FILE AND PROCESS DO INVERSION WITH N CYCLES 


	if(configCrontrolFile.NumberOfCycles<0){
		// read fits or per 
      
		AllocateMemoryDerivedSynthesis(nlambda);
		if(strcmp(file_ext(configCrontrolFile.ObservedProfiles),PER_FILE)==0){ // invert only per file
			float * spectroPER = calloc(nlambda*NPARMS,sizeof(float));
			FILE * fReadSpectro;
			char * line = NULL;
			size_t len = 0;
			ssize_t read;
			fReadSpectro = fopen(configCrontrolFile.ObservedProfiles, "r");
			
			int contLine=0;
			if (fReadSpectro == NULL)
			{
				printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
				printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", configCrontrolFile.ObservedProfiles);
				fclose(fReadSpectro);
				exit(EXIT_FAILURE);
			}
			float aux1, aux2,aux3,aux4,aux5,aux6;
			while ((read = getline(&line, &len, fReadSpectro)) != -1 && contLine<nlambda) {
				//sscanf(line,"%e %e %e %e %e %e",&indexLine,&dummy,&spectroPER[contLine], &spectroPER[contLine + nlambda], &spectroPER[contLine + nlambda * 2], &spectroPER[contLine + nlambda * 3]);
				sscanf(line,"%e %e %e %e %e %e",&aux1,&aux2,&aux3,&aux4,&aux5,&aux6);
				spectroPER[contLine] = aux3;
				spectroPER[contLine + nlambda] = aux4;
				spectroPER[contLine + nlambda * 2] = aux5;
				spectroPER[contLine + nlambda * 3] = aux6;
				contLine++;
			}
			fclose(fReadSpectro);

      	

			Init_Model initModel;
			initModel.eta0 = 0;
			initModel.mac = 0;
			initModel.dopp = 0;
			initModel.aa = 0;
			initModel.alfa = 0; //0.38; //stray light factor
			initModel.S1 = 0;
			//invert with classical estimates
			estimacionesClasicas(wlines[1], vLambda, nlambda, spectroPER, &initModel,0);
			// save model to file
			char nameAuxOutputModel [4096];
			if(configCrontrolFile.ObservedProfiles[0]!='\0')
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.ObservedProfiles));
			else
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.InitialGuessModel));
				

			strcat(nameAuxOutputModel,"_model_ce");
			strcat(nameAuxOutputModel,MOD_FILE);
			FILE * fptr = fopen(nameAuxOutputModel, "w");
			if(fptr!=NULL){
				fprintf(fptr,"MODEL_ETHA0: %lf",initModel.eta0);
				fprintf(fptr,"\nMODEL_B: %lf",initModel.B);
				fprintf(fptr,"\nMODEL_VLOS: %lf",initModel.vlos);
				fprintf(fptr,"\nMODEL_LAMBDADOPP: %lf",initModel.dopp);
				fprintf(fptr,"\nMODEL_AA: %lf",initModel.aa);
				fprintf(fptr,"\nMODEL_GM: %lf",initModel.gm);
				fprintf(fptr,"\nMODEL_AZI: %lf",initModel.az);
				fprintf(fptr,"\nMODEL_S0: %lf",initModel.S0);
				fprintf(fptr,"\nMODEL_S1: %lf",initModel.S1);
				fprintf(fptr,"\nMODEL_MAC: %lf",initModel.mac);
				fprintf(fptr,"\nMODEL_ALFA: %lf",initModel.alfa);
				fprintf(fptr,"\nNUMBER ITERATIONS TO CONVERGE: %d",0);
				fprintf(fptr,"\nChisqrf %le",0.0);
				fprintf(fptr,"\n\n");
				fclose(fptr);
			}
			else{
				printf("\n ¡¡¡¡¡ ERROR: OUTPUT MODEL FILE CAN NOT BE OPENED: %s \n !!!!!",nameAuxOutputModel);
			}			
			free(spectroPER);
		}
		else if(strcmp(file_ext(configCrontrolFile.ObservedProfiles),FITS_FILE)==0){ // invert image from fits file 
			fitsImage = readFitsSpectroImage(configCrontrolFile.ObservedProfiles,0);
			//fitsImage = readFitsSpectroImageRectangular(configCrontrolFile.ObservedProfiles,&configCrontrolFile,0);
			// ALLOCATE MEMORY FOR STORE THE RESULTS 
			int indexPixel = 0;
			vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
			vChisqrf = calloc (fitsImage->numPixels , sizeof(float));
			vNumIter = calloc (fitsImage->numPixels , sizeof(int));
			for(indexPixel = 0; indexPixel < fitsImage->numPixels; indexPixel++){

				//Initial Model
				Init_Model initModel;
				initModel.eta0 = INITIAL_MODEL.eta0;
				initModel.B = INITIAL_MODEL.B; //200 700
				initModel.gm = INITIAL_MODEL.gm;
				initModel.az = INITIAL_MODEL.az;
				initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
				initModel.mac = INITIAL_MODEL.mac;
				initModel.dopp = INITIAL_MODEL.dopp;
				initModel.aa = INITIAL_MODEL.aa;
				initModel.alfa = INITIAL_MODEL.alfa; //0.38; //stray light factor
				initModel.S0 = INITIAL_MODEL.S0;
				initModel.S1 = INITIAL_MODEL.S1;
				estimacionesClasicas(wlines[1],vLambda, nlambda, fitsImage->pixels[indexPixel].spectro, &initModel,0);
				vModels[indexPixel] = initModel;

			}
			char nameAuxOutputModel [4096];
			if(configCrontrolFile.ObservedProfiles[0]!='\0')
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.ObservedProfiles));
			else
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.InitialGuessModel));
			strcat(nameAuxOutputModel,"_model_ce");
			strcat(nameAuxOutputModel,FITS_FILE);			
			if(!writeFitsImageModels(nameAuxOutputModel,fitsImage->rows,fitsImage->cols,vModels,vChisqrf,vNumIter,configCrontrolFile.saveChisqr)){
					printf("\n ERROR WRITING FILE OF MODELS: %s",nameAuxOutputModel);
					//slog_error(0,"\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
			}
			free(vModels);
			free(vChisqrf);
			free(vNumIter);
		}
		else{
			printf("\n OBSERVED PROFILES DOESN'T HAVE CORRECT EXTENSION  .PER or .FITS ");
			exit(EXIT_FAILURE);
		}
	}
	else if(configCrontrolFile.NumberOfCycles==0){ // synthesis
	
      if(access(configCrontrolFile.StrayLightFile,F_OK)!=-1){ //  IF NOT EMPTY READ stray light file 
         slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,nlambda);
         printf("\n STRAY LIGHT READ \n");
      }
      else{
         printf("\n STRAY LIGHT NOT USED \n");
      }      
      Init_Model initModel;
      initModel.eta0 = INITIAL_MODEL.eta0;
      initModel.B = INITIAL_MODEL.B; //200 700
      initModel.gm = INITIAL_MODEL.gm;
      initModel.az = INITIAL_MODEL.az;
      initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
      initModel.mac = INITIAL_MODEL.mac;
      initModel.dopp = INITIAL_MODEL.dopp;
      initModel.aa = INITIAL_MODEL.aa;
      initModel.alfa = INITIAL_MODEL.alfa; //0.38; //stray light factor
      initModel.S0 = INITIAL_MODEL.S0;
      initModel.S1 = INITIAL_MODEL.S1;
      printf("\n MODEL ATMOSPHERE: \n");
      printf("\n ETA0: %lf",initModel.eta0);
      printf("\n B: %lf",initModel.B);
      printf("\n vlos: %lf",initModel.vlos);
      printf("\n dopp: %lf",initModel.dopp);
      printf("\n aa: %lf",initModel.aa);
      printf("\n gm: %lf",initModel.gm);
      printf("\n az: %lf",initModel.az);
      printf("\n S0: %lf",initModel.S0);
      printf("\n S1: %lf",initModel.S1);      
      printf("\n mac: %lf",initModel.mac);
      printf("\n alfa: %lf",initModel.alfa);
		printf("\n");    

      
      
      AllocateMemoryDerivedSynthesis(nlambda);

		// synthesis
      mil_sinrf(cuantic, &initModel, wlines, vLambda, nlambda, spectra, AH, slight,spectra_mac, configCrontrolFile.ConvolveWithPSF);
      me_der(cuantic, &initModel, wlines, vLambda, nlambda, d_spectra, spectra_mac, spectra, AH, slight, configCrontrolFile.ConvolveWithPSF);	

		// in this case basenamefile is from initmodel
		char nameAux [4096];
		if(configCrontrolFile.ObservedProfiles[0]!='\0')
			strcpy(nameAux,get_basefilename(configCrontrolFile.ObservedProfiles));
		else
			strcpy(nameAux,get_basefilename(configCrontrolFile.InitialGuessModel));		
		strcat(nameAux,PER_FILE);
		FILE *fptr = fopen(nameAux, "w");
		if(fptr!=NULL){
			int kk;
			for (kk = 0; kk < nlambda; kk++)
			{
				//fprintf(fptr,"%d\t%f\t%e\t%e\t%e\t%e\n", indexLine, (vLambda[kk]-configCrontrolFile.CentralWaveLenght)*1000, spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
				fprintf(fptr,"%d\t%f\t%e\t%e\t%e\t%e\n", indexLine, vLambda[kk], spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
			}
			fclose(fptr);
			printf("\n*******************************************************************************************");
			printf("\n******************SYNTHESIS DONE: %s",nameAux);
			printf("\n*******************************************************************************************\n");
		}
		else{
			printf("\n ERROR !!! The output file can not be open: %s",nameAux);
		}

		/*int number_parametros = 0;
		for (number_parametros = 0; number_parametros < NTERMS; number_parametros++)
		{
			strcpy(nameAux,get_basefilename(configCrontrolFile.InitialGuessModel));
			strcat(nameAux,"_C_");
			char extension[10];
			sprintf(extension, "%d%s", number_parametros,".per");
			strcat(nameAux,extension);
			FILE *fptr = fopen(nameAux, "w");
			//printf("\n FUNCION RESPUESTA: %d \n",number_parametros);
			int kk;
			for (kk = 0; kk < nlambda; kk++)
			{
			fprintf(fptr,"1\t%lf\t%le\t%le\t%le\t%le\n", vLambda[kk],
			d_spectra[kk + nlambda * number_parametros],
			d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS],
			d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 2],
			d_spectra[kk + nlambda * number_parametros + nlambda * NTERMS * 3]);
			}
			fclose(fptr);
		}
		printf("\n");*/

	}
	else{ // INVERT PIXEL FROM PER FILE OR IMAGE FROM FITS FILE 
		if(strcmp(file_ext(configCrontrolFile.ObservedProfiles),PER_FILE)==0){ // invert only per file
			float * spectroPER = calloc(nlambda*NPARMS,sizeof(float));
			FILE * fReadSpectro;
			char * line = NULL;
			size_t len = 0;
			ssize_t read;
			fReadSpectro = fopen(configCrontrolFile.ObservedProfiles, "r");
			
			int contLine=0;
			if (fReadSpectro == NULL)
			{
				printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
				printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", configCrontrolFile.ObservedProfiles);
				fclose(fReadSpectro);
				exit(EXIT_FAILURE);
			}
			
			float aux1, aux2, aux3, aux4, aux5, aux6;
			while ((read = getline(&line, &len, fReadSpectro)) != -1 && contLine<nlambda) {
				sscanf(line,"%e %e %e %e %e %e",&aux1,&aux2,&aux3,&aux4,&aux5,&aux6);
				spectroPER[contLine] = aux3;
				spectroPER[contLine + nlambda] = aux4;
				spectroPER[contLine + nlambda * 2] = aux5;
				spectroPER[contLine + nlambda * 3] = aux6;
				contLine++;
			}
			fclose(fReadSpectro);

      
      	
			AllocateMemoryDerivedSynthesis(nlambda);
			Init_Model initModel;
			initModel.eta0 = INITIAL_MODEL.eta0;
			initModel.B = INITIAL_MODEL.B; //200 700
			initModel.gm = INITIAL_MODEL.gm;
			initModel.az = INITIAL_MODEL.az;
			initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
			initModel.mac = INITIAL_MODEL.mac;
			initModel.dopp = INITIAL_MODEL.dopp;
			initModel.aa = INITIAL_MODEL.aa;
			initModel.alfa = INITIAL_MODEL.alfa; //0.38; //stray light factor
			initModel.S0 = INITIAL_MODEL.S0;
			initModel.S1 = INITIAL_MODEL.S1;

      	int numIter;
      	lm_mils(cuantic, wlines, vLambda, nlambda, spectroPER, nlambda, &initModel, spectra, &chisqrf, slight, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
               configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&numIter);

			// SAVE OUTPUT MODEL 
			char nameAuxOutputModel [4096];
			/*strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.ObservedProfiles));
			strcat(nameAuxOutputModel,OUTPUT_MOD_TXT_EXT);*/
			if(configCrontrolFile.ObservedProfiles[0]!='\0')
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.ObservedProfiles));
			else
				strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.InitialGuessModel));
							
			strcat(nameAuxOutputModel,"_model");
			strcat(nameAuxOutputModel,MOD_FILE);

			FILE *fptr = fopen(nameAuxOutputModel, "w");
			if(fptr!=NULL){
				fprintf(fptr,"MODEL_ETHA0: %lf\n",initModel.eta0);
				fprintf(fptr,"MODEL_B: %lf\n",initModel.B);
				fprintf(fptr,"MODEL_VLOS: %lf\n",initModel.vlos);
				fprintf(fptr,"MODEL_LAMBDADOPP: %lf\n",initModel.dopp);
				fprintf(fptr,"MODEL_AA: %lf\n",initModel.aa);
				fprintf(fptr,"MODEL_GM: %lf\n",initModel.gm);
				fprintf(fptr,"MODEL_AZI: %lf\n",initModel.az);
				fprintf(fptr,"MODEL_S0: %lf\n",initModel.S0);
				fprintf(fptr,"MODEL_S1: %lf\n",initModel.S1);
				fprintf(fptr,"MODEL_MAC: %lf\n",initModel.mac);
				fprintf(fptr,"MODEL_ALFA: %lf\n",initModel.alfa);
				fprintf(fptr,"NUMBER ITERATIONS TO CONVERGE: %d\n",numIter);
				fprintf(fptr,"Chisqrf %le\n",chisqrf);
				fprintf(fptr,"\n\n");
				fclose(fptr);
				printf("\n*******************************************************************************************");
				printf("\n******************INVERTED MODEL SAVED IN FILE: %s",nameAuxOutputModel);
				printf("\n*******************************************************************************************\n");				
			}
			else{
				printf("\n ¡¡¡¡¡ ERROR: OUTPUT MODEL FILE CAN NOT BE OPENED\n !!!!! ");
			}


			// SAVE OUTPUT ADJUST SYNTHESIS PROFILES 
			if(configCrontrolFile.SaveSynthesisAdjusted){
				char nameAuxOutputStokes [4096];
				//strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.ObservedProfiles));
				//strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.InitialGuessModel));
				if(configCrontrolFile.ObservedProfiles[0]!='\0')
					strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.ObservedProfiles));
				else
					strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.InitialGuessModel));				
				strcat(nameAuxOutputStokes,STOKES_PER_EXT);
				FILE *fptr = fopen(nameAuxOutputStokes, "w");
				if(fptr!=NULL){
			      //printf("\n valores de spectro sintetizado\n");
					int kk;
					for (kk = 0; kk < nlambda; kk++)
					{
						fprintf(fptr,"%d\t%le\t%le\t%le\t%le\t%le\n", indexLine, vLambda[kk]-configCrontrolFile.CentralWaveLenght, spectra[kk], spectra[kk + nlambda], spectra[kk + nlambda * 2], spectra[kk + nlambda * 3]);
					}
					//printf("\nVALORES DE LAS FUNCIONES RESPUESTA \n");
					fclose(fptr);
					printf("\n*******************************************************************************************");
					printf("\n******************SPECTRUM SYNTHESIS ADJUSTED SAVED IN FILE: %s",nameAuxOutputStokes);
					printf("\n*******************************************************************************************\n\n");					
				}
				else{
					printf("\n ¡¡¡¡¡ ERROR: OUTPUT SYNTHESIS PROFILE ADJUSTED FILE CAN NOT BE OPENED\n !!!!! ");
				}
			}
			free(spectroPER);	
		}
		else if(strcmp(file_ext(configCrontrolFile.ObservedProfiles),FITS_FILE)==0){ // invert image from fits file 

			// READ PIXELS FROM IMAGE 
			PRECISION timeReadImage;
			clock_t t;
			t = clock();
			
			
			//fitsImage = readFitsSpectroImageRectangular(configCrontrolFile.ObservedProfiles,&configCrontrolFile,0);
			fitsImage = readFitsSpectroImage(nameInputFileSpectra,0);
			t = clock() - t;
			timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			
			printf("\n\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
			//slog_info(0,"\n\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage);

			if(fitsImage!=NULL){
				FitsImage * imageStokesAdjust = NULL;
				if(configCrontrolFile.SaveSynthesisAdjusted){
					imageStokesAdjust = malloc(sizeof(FitsImage));
					imageStokesAdjust->rows = fitsImage->rows;
					imageStokesAdjust->cols = fitsImage->cols;
					imageStokesAdjust->nLambdas = fitsImage->nLambdas;
					imageStokesAdjust->numStokes = fitsImage->numStokes;
					imageStokesAdjust->pos_col = fitsImage->pos_col;
					imageStokesAdjust->pos_row = fitsImage->pos_row;
					imageStokesAdjust->pos_lambda = fitsImage->pos_lambda;
					imageStokesAdjust->pos_stokes_parameters = fitsImage->pos_stokes_parameters;
					imageStokesAdjust->numPixels = fitsImage->numPixels;
					imageStokesAdjust->pixels = calloc(imageStokesAdjust->numPixels, sizeof(vpixels));
					for( i=0;i<imageStokesAdjust->numPixels;i++){
						imageStokesAdjust->pixels[i].spectro = calloc ((imageStokesAdjust->numStokes*imageStokesAdjust->nLambdas),sizeof(float));
					}
				}				
				// check if read stray light
				if(access(configCrontrolFile.StrayLightFile,F_OK)!=-1){ //  IF NOT EMPTY READ stray light file 
					slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,nlambda);
				}

				//***************************************** INIT MEMORY WITH SIZE OF LAMBDA ****************************************************//
				AllocateMemoryDerivedSynthesis(nlambda);
				int indexPixel = 0;

				// ALLOCATE MEMORY FOR STORE THE RESULTS 

				vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
				vChisqrf = calloc (fitsImage->numPixels , sizeof(float));
				vNumIter = calloc (fitsImage->numPixels , sizeof(int));
				t = clock();
				printf("\n***********************  PROGRESS INVERSION *******************************\n\n");
				//slog_info(0,"\n***********************  PROGRESS INVERSION *******************************\n\n");

				for(indexPixel = 0; indexPixel < fitsImage->numPixels; indexPixel++){

					//Initial Model
					Init_Model initModel;
					initModel.eta0 = INITIAL_MODEL.eta0;
					initModel.B = INITIAL_MODEL.B; //200 700
					initModel.gm = INITIAL_MODEL.gm;
					initModel.az = INITIAL_MODEL.az;
					initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
					initModel.mac = INITIAL_MODEL.mac;
					initModel.dopp = INITIAL_MODEL.dopp;
					initModel.aa = INITIAL_MODEL.aa;
					initModel.alfa = INITIAL_MODEL.alfa; //0.38; //stray light factor
					initModel.S0 = INITIAL_MODEL.S0;
					initModel.S1 = INITIAL_MODEL.S1;
					
					// CLASSICAL ESTIMATES TO GET B, GAMMA
					estimacionesClasicas(wlines[1], vLambda, nlambda, fitsImage->pixels[indexPixel].spectro, &initModel,1);
					if (isnan(initModel.B))
						initModel.B = 1;
					if (isnan(initModel.vlos))
						initModel.vlos = 1e-3;
					if (isnan(initModel.gm))
						initModel.gm = 1;						
					if (isnan(initModel.az))
						initModel.az = 1;
					// INVERSION RTE

					PRECISION * slightPixel;
					if(slight==NULL) 
						slightPixel = NULL;
					else{
						if(dimStrayLight==nlambda) 
							slightPixel = slight;
						else 
							slightPixel = slight+nlambda*indexPixel;
					}
					lm_mils(cuantic, wlines, vLambda, nlambda, fitsImage->pixels[indexPixel].spectro, nlambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
							configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);						
					
					vModels[indexPixel] = initModel;
					if(configCrontrolFile.SaveSynthesisAdjusted){
						int kk;
						for (kk = 0; kk < (nlambda * NPARMS); kk++)
						{
							imageStokesAdjust->pixels[indexPixel].spectro[kk] = spectra[kk] ;
						}						
					}					
					//printf ("\t\t %.2f seconds -- %.2f %%\r",  ((PRECISION)(clock() - t)/CLOCKS_PER_SEC) , ((indexPixel*100.)/fitsImage->numPixels));
				}
				t = clock() - t;
				timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n FINISH EXECUTION OF INVERSION: %f seconds to execute \n", timeReadImage);
				//slog_info(0,"\n FINISH EXECUTION OF INVERSION: %f seconds to execute \n", timeReadImage);
				
				char nameAuxOutputModel [4096];
				if(configCrontrolFile.ObservedProfiles[0]!='\0')
					strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.ObservedProfiles));
				else
					strcpy(nameAuxOutputModel,get_basefilename(configCrontrolFile.InitialGuessModel));				

				strcat(nameAuxOutputModel,MOD_FITS);
				if(!writeFitsImageModels(nameAuxOutputModel,fitsImage->rows,fitsImage->cols,vModels,vChisqrf,vNumIter,configCrontrolFile.saveChisqr)){
						printf("\n ERROR WRITING FILE OF MODELS: %s",nameAuxOutputModel);
						//slog_error(0,"\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
				}
				// PROCESS FILE OF SYNTETIC PROFILES

				if(configCrontrolFile.SaveSynthesisAdjusted){
					// WRITE SINTHETIC PROFILES TO FITS FILE
					char nameAuxOutputStokes [4096];
					//strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.ObservedProfiles));
					//strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.InitialGuessModel));
					if(configCrontrolFile.ObservedProfiles[0]!='\0')
						strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.ObservedProfiles));
					else
						strcpy(nameAuxOutputStokes,get_basefilename(configCrontrolFile.InitialGuessModel));					
					strcat(nameAuxOutputStokes,STOKES_FIT_EXT);
					if(!writeFitsImageProfiles(nameAuxOutputStokes,nameInputFileSpectra,imageStokesAdjust)){
						printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",nameOutputFilePerfiles);
					}
				}
				if(configCrontrolFile.SaveSynthesisAdjusted)
					free(imageStokesAdjust);
				free(vModels);
				free(vChisqrf);
				free(vNumIter);					
			}
			else{
				printf("\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
			}			

			printf(" \n***********************  IMAGE INVERSION DONE, CLEANING MEMORY *********************\n");
			freeFitsImage(fitsImage);
		}
		else{
			printf("\n OBSERVED PROFILES DOESN'T HAVE CORRECT EXTENSION  .PER or .FITS ");
			exit(EXIT_FAILURE);
		}
	}

	fftw_free(inFilterMAC);
	fftw_free(outFilterMAC);
	fftw_destroy_plan(planFilterMAC);
	fftw_free(inFilterMAC_DERIV);
	fftw_free(outFilterMAC_DERIV);
	fftw_destroy_plan(planFilterMAC_DERIV);
	fftw_free(inSpectraFwMAC);
	fftw_free(outSpectraFwMAC);
	fftw_destroy_plan(planForwardMAC);
	fftw_free(inSpectraBwMAC);
	fftw_free(outSpectraBwMAC);
	fftw_destroy_plan(planBackwardMAC);

	if(configCrontrolFile.ConvolveWithPSF){
		fftw_free(inSpectraFwPSF);
		fftw_free(outSpectraFwPSF);
		fftw_destroy_plan(planForwardPSF);
		fftw_free(inSpectraBwPSF);
		fftw_free(outSpectraBwPSF);
		fftw_destroy_plan(planBackwardPSF);

		fftw_free(fftw_G_PSF);
		fftw_free(fftw_G_MAC_PSF);
		fftw_free(fftw_G_MAC_DERIV_PSF);

		fftw_free(inPSF_MAC);
		fftw_free(inMulMacPSF);
		fftw_free(inPSF_MAC_DERIV);
		fftw_free(inMulMacPSFDeriv);
		fftw_free(outConvFilters);
		fftw_free(outConvFiltersDeriv);	

		fftw_destroy_plan(planForwardPSF_MAC);
		fftw_destroy_plan(planForwardPSF_MAC_DERIV);
		fftw_destroy_plan(planBackwardPSF_MAC);
		fftw_destroy_plan(planBackwardPSF_MAC_DERIV);		
	}

	free(cuantic);
	free(wlines);
	FreeMemoryDerivedSynthesis();
	if(G!=NULL) free(G);
	gsl_eigen_symmv_free (workspace);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	return 0;
}
