
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
// CMILOS v1.0 (2019)
// RTE INVERSION C - MPI code for SOPHI (based on the ILD code MILOS by D. Orozco)
// Juanp && Manu (IAA-CSIC)


#include "mpi.h"
#include <time.h>
#include "defines.h"
//#include "nrutil.h"
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "readConfig.h"
#include "lib.h"
#include "milosUtils.h"
#include <unistd.h>
#include <complex.h>
#include <fftw3.h> //always after complex.h
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>



// ***************************** FUNCTIONS TO READ FITS FILE *********************************************************


long long int c1, c2, cd, semi, c1a, c2a, cda; //variables of 64 bits to read whatch cycles 
long long int c1total, c2total, cdtotal, ctritotal;

Cuantic *cuantic; // Global variable with cuantic information 


PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

PRECISION *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
PRECISION *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;
PRECISION *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
PRECISION *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
PRECISION *dfi, *dshi;
PRECISION CC, CC_2, sin_gm, azi_2, sinis, cosis, cosis_2, cosi, sina, cosa, sinda, cosda, sindi, cosdi, sinis_cosa, sinis_sina;
PRECISION *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
PRECISION *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
PRECISION *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
PRECISION *parcial1, *parcial2, *parcial3;
PRECISION *nubB, *nupB, *nurB;
PRECISION **uuGlobalInicial;
PRECISION **HGlobalInicial;
PRECISION **FGlobalInicial;
PRECISION *perfil_instrumental;
PRECISION *G,*GMAC;
PRECISION *interpolatedPSF;
PRECISION AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
PRECISION * opa;
int FGlobal, HGlobal, uuGlobal;

PRECISION *d_spectra, *spectra, *spectra_mac;


// GLOBAL variables to use for FFT calculation 

fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
fftw_complex * inSpectraFwMAC, *inSpectraBwMAC, *outSpectraFwMAC, *outSpectraBwMAC;
fftw_plan planForwardPSF, planBackwardPSF;
fftw_plan planForwardMAC, planBackwardMAC;
fftw_complex * inFilterMAC, * inFilterMAC_DERIV, * outFilterMAC, * outFilterMAC_DERIV;
fftw_plan planFilterMAC, planFilterMAC_DERIV;
fftw_complex * fftw_G_PSF;

fftw_complex * fftw_G_PSF, * fftw_G_MAC_PSF, * fftw_G_MAC_DERIV_PSF;
fftw_complex * inPSF_MAC, * inMulMacPSF, * inPSF_MAC_DERIV, *inMulMacPSFDeriv, *outConvFilters, * outConvFiltersDeriv;
fftw_plan planForwardPSF_MAC, planForwardPSF_MAC_DERIV,planBackwardPSF_MAC, planBackwardPSF_MAC_DERIV;

//Convolutions values
int sizeG = 0;
PRECISION FWHM = 0;

ConfigControl configCrontrolFile;

_Complex PRECISION *z,* zden, * zdiv;

int main(int argc, char **argv)
{
	int i, j;  // indexes 
	
	int indexLine; // index to identify central line to read it 
	PRECISION initialLambda, step, finalLambda;
	// INIT MPI  PROGRAM 
	int numProcs, idProc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &idProc);

	/** DEFINE MPI TYPE TO SEND MODELS **/
	MPI_Datatype mpiInitModel;
	const int nitemsStructInitModel = 11;
	int blocklenghtInitModel [11] = {1,1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype typesInitModel [11] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};					
	MPI_Aint offsetsInitModel [11];
	offsetsInitModel[0] = offsetof(Init_Model, eta0);
	offsetsInitModel[1] = offsetof(Init_Model, B);
	offsetsInitModel[2] = offsetof(Init_Model, vlos);
	offsetsInitModel[3] = offsetof(Init_Model, dopp);
	offsetsInitModel[4] = offsetof(Init_Model, aa);
	offsetsInitModel[5] = offsetof(Init_Model, gm);
	offsetsInitModel[6] = offsetof(Init_Model, az);
	offsetsInitModel[7] = offsetof(Init_Model, S0);
	offsetsInitModel[8] = offsetof(Init_Model, S1);
	offsetsInitModel[9] = offsetof(Init_Model, mac);
	offsetsInitModel[10] = offsetof(Init_Model, alfa);
	MPI_Type_create_struct(nitemsStructInitModel, blocklenghtInitModel, offsetsInitModel, typesInitModel, &mpiInitModel);
	MPI_Type_commit(&mpiInitModel);


	MPI_Datatype mpiName;
	const int nItemsStructName = 1;
	int blocklenghName [1] = {PATH_MAX};
	MPI_Datatype typesName [1] = {MPI_CHAR};
	MPI_Aint offsetName [1];
	offsetName[0] = offsetof(nameFile,name) ;
	MPI_Type_create_struct(nItemsStructName,blocklenghName,offsetName,typesName,&mpiName);
	MPI_Type_commit(&mpiName);

	MPI_Request mpiRequestSpectro, mpiRequestLambda,mpiRequestName;

	/************************************/
	const int root=0;	
	
	// FINISH STARTING PROGRAM 

	PRECISION *wlines;
	int nlambda, numPixels, indexPixel;
	
	int posCENTRAL_WL; // position central wl in file of LINES
	Init_Model INITIAL_MODEL;
	PRECISION * deltaLambda, * PSF;
	int N_SAMPLES_PSF;
	
	
	// CONFIGURACION DE PARAMETROS A INVERTIR
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	//----------------------------------------------

	PRECISION * slight = NULL;
	
	int dimStrayLight;
	int numberOfFileSpectra;

   nameFile * vInputFileSpectra;
	nameFile * vInputFileSpectraParalell = NULL;
	nameFile * vOutputNameModels;
	nameFile * vOutputNameModelsParalell = NULL;
	nameFile * vOutputNameSynthesisAdjusted;
	nameFile * vOutputNameSynthesisAdjustedParallel = NULL;
	nameFile * vInputFileSpectraLocal;
	nameFile * vOutputNameModelsLocal;
	nameFile * vOutputNameSynthesisAdjustedLocal;

	char * baseNameOutputObserved;
	const char	* nameInputFilePSF ;

   FitsImage * fitsImage = NULL;
	PRECISION  dat[7];

	double local_start, local_finish, local_elapsed, elapsed;
	double local_start_execution, local_finish_execution, local_elapsed_execution, elapsed_execution;
	double local_start_scatter, local_finish_scatter, local_elapsed_scatter, elapsed_scatter;
	double local_start_gather, local_finish_gather, local_elapsed_gather, elapsed_gather;
	
	Init_Model * resultsInitModel;
	Init_Model * resultsInitModelTotal;
	PRECISION * chisqrfTotal, *vChisqrf, chisqrf;
	int * vNumIter, * vNumIterTotal; // to store the number of iterations used to converge for each pixel
	float  * vSpectraSplit, * vLambdaSplit, * vSpectraAdjustedSplit, * vSpectraAjustedTotal;

	int sendcountsPixels [numProcs] ; // array describing how many elements to send to each process
	int sendcountsSpectro [numProcs];
	int sendcountsLambda [numProcs];
	int sendcountsNameInputFiles [numProcs];  // how many files per process
	int displsPixels [numProcs];  // array describing the displacements where each segment begins
	int displsSpectro [numProcs];
	int displsLambda [numProcs];
	int displsNameInputFiles [numProcs]; // how many 

	/********************* Read data input from file ******************************/

	loadInitialValues(&configCrontrolFile);
	//readParametersFileInput(argv[1],&configCrontrolFile,(idProc==root));
	if(!readInitFile(argv[1],&configCrontrolFile,(idProc==root))){
		if(idProc==root)
			printf("\n\n ¡¡¡ ERROR READING INIT FILE !!! \n\n");
		exit(EXIT_FAILURE);
	}
	
	// check if type of file is FITS, in other case exit 
	if(strcmp(configCrontrolFile.typeInputStokes,"fits")!=0){
		if(idProc==root)
			printf("\n ERROR, the files in parallel version must be in FITS file only.\n");
		exit(EXIT_FAILURE);
	}

	// CHECK IF ONLY RECEIVED ONE FILE AND THE EXTENSION IS .FITS 
	if(configCrontrolFile.t1 == 0 && configCrontrolFile.t2 ==0){ // then process only one file
		if(strcmp(file_ext(configCrontrolFile.ObservedProfiles),FITS_FILE)!=0){ 
			if(idProc==root){
				printf("\n**********************************************************\n");
				printf("\nERROR, without specify timeseries the value of control parameter 'Observed Profiles' must be a fits file\n");
				printf("\n**********************************************************\n");
			}
			exit(EXIT_FAILURE);
		}
	}
	
	nameInputFilePSF = configCrontrolFile.PSFFile;
	FWHM = configCrontrolFile.FWHM;

	/***************** READ INIT MODEL ********************************/
	if(!readInitialModel(&INITIAL_MODEL,configCrontrolFile.InitialGuessModel)){
		printf("\n\n ¡¡¡ ERROR READING INIT MODEL !!! \n\n");
		exit(EXIT_FAILURE);
	}

	/***************** READ WAVELENGHT FROM GRID OR FITS ********************************/
	PRECISION * vGlobalLambda;

	if(configCrontrolFile.useMallaGrid){ // read lambda from grid file
		indexLine = readMallaGrid(configCrontrolFile.MallaGrid, &initialLambda, &step, &finalLambda, (idProc==root));      
		nlambda = ((finalLambda-initialLambda)/step)+1;
		// pass to armstrong 
		initialLambda = initialLambda/1000;
		step = step/1000;
		finalLambda = finalLambda/1000;
		vGlobalLambda = calloc(nlambda,sizeof(PRECISION));
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(configCrontrolFile.AtomicParametersFile,dat,indexLine,(idProc==root));
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n CUANTIC LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVELENGHT: %f",configCrontrolFile.CentralWaveLenght);
			exit(1);
		}
		vGlobalLambda[0]=configCrontrolFile.CentralWaveLenght+(initialLambda);
   		for(i=1;i<nlambda;i++){
        	vGlobalLambda[i]=vGlobalLambda[i-1]+step;
     	}
	}
	else{
		vGlobalLambda = readFitsLambdaToArray(configCrontrolFile.WavelengthFile,0,0,&indexLine,&nlambda);
		if(vGlobalLambda==NULL){
			printf("\n FILE WITH WAVELENGHT HAS NOT BEEN READ PROPERLY, please check it.\n");
			free(vGlobalLambda);
			exit(EXIT_FAILURE);
		}
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(configCrontrolFile.AtomicParametersFile,dat,indexLine,(idProc==root));
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n CUANTIC LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
			exit(1);
		}
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	/*********************************************** INITIALIZE VARIABLES  *********************************/

	CC = PI / 180.0;
	CC_2 = CC * 2;
	
	wlines = (PRECISION *)calloc(2, sizeof(PRECISION));
	wlines[0] = 1;
	wlines[1] = configCrontrolFile.CentralWaveLenght;

	numPixels=0;	
	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/
	MPI_Barrier(MPI_COMM_WORLD);
	cuantic = create_cuantic(dat,(idProc==root));
	InitializePointerShareCalculation();

	/**************************************** READ FITS  STRAY LIGHT ******************************/
	MPI_Barrier(MPI_COMM_WORLD);
	if(access(configCrontrolFile.StrayLightFile,F_OK)!=-1){ //  IF NOT EMPTY READ stray light file 
		slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,nlambda);
	}

	
	// ************************** DEFINE PLANS TO EXECUTE MACROTURBULENCE IF NECESSARY **********************************************//
	// MACROTURBULENCE PLANS
	
	inFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC = fftw_plan_dft_1d(nlambda, inFilterMAC, outFilterMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
	inFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC_DERIV = fftw_plan_dft_1d(nlambda, inFilterMAC_DERIV, outFilterMAC_DERIV, FFT_FORWARD, FFTW_EXHAUSTIVE);


	inSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planForwardMAC = fftw_plan_dft_1d(nlambda, inSpectraFwMAC, outSpectraFwMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
	inSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
	planBackwardMAC = fftw_plan_dft_1d(nlambda, inSpectraBwMAC, outSpectraBwMAC, FFT_BACKWARD, FFTW_EXHAUSTIVE);

	// ********************************************* IF PSF HAS BEEN SELECTEC IN TROL READ PSF FILE OR CREATE GAUSSIAN FILTER ***********//
	if(configCrontrolFile.ConvolveWithPSF){

		if(configCrontrolFile.FWHM > 0){
			G = fgauss_WL(FWHM,vGlobalLambda[1]-vGlobalLambda[0],vGlobalLambda[0],vGlobalLambda[nlambda/2],nlambda,&sizeG);
		}
		else if(access(nameInputFilePSF,F_OK) != -1){
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
					readPSFFile(deltaLambda,PSF,nameInputFilePSF);
					PRECISION * fInterpolated = calloc(nlambda,sizeof(PRECISION));
					interpolationLinearPSF(deltaLambda,  PSF, vGlobalLambda ,configCrontrolFile.CentralWaveLenght, N_SAMPLES_PSF,fInterpolated, nlambda);						
				}
				else{
					//PRECISION * fgauss_WL(PRECISION FWHM, PRECISION step_between_lw, PRECISION lambda0, PRECISION lambdaCentral, int nLambda, int * sizeG)
					G = fgauss_WL(FWHM,vGlobalLambda[1]-vGlobalLambda[0],vGlobalLambda[0],vGlobalLambda[nlambda/2],nlambda,&sizeG);
				}
		}

		
		//PSF FILTER PLANS 
		inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF = fftw_plan_dft_1d(nlambda, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
		planBackwardPSF = fftw_plan_dft_1d(nlambda, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_EXHAUSTIVE);

		fftw_complex * in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		int i;
		for (i = 0; i < nlambda; i++)
		{
			in[i] = G[i] + 0 * _Complex_I;
		}
		fftw_G_PSF = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_plan p = fftw_plan_dft_1d(nlambda, in, fftw_G_PSF, FFT_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		for (i = 0; i < nlambda; i++)
		{
			fftw_G_PSF[i] = fftw_G_PSF[i] / nlambda;
		}
		fftw_destroy_plan(p);
		fftw_free(in);
		
		inPSF_MAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC = fftw_plan_dft_1d(nlambda, inPSF_MAC, fftw_G_MAC_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inMulMacPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFilters = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC = fftw_plan_dft_1d(nlambda, inMulMacPSF, outConvFilters, FFT_BACKWARD, FFTW_EXHAUSTIVE);


		inPSF_MAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_DERIV_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inPSF_MAC_DERIV, fftw_G_MAC_DERIV_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inMulMacPSFDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFiltersDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inMulMacPSFDeriv, outConvFiltersDeriv, FFT_BACKWARD, FFTW_EXHAUSTIVE);			

	}

	/*****************************************************************************************************/
	MPI_Barrier(MPI_COMM_WORLD);
	// ROOT PROCESS READ IMAGE FROM FILE TO KNOW LIST OF FILES
	if(idProc==root){

      	// CHECK IF INPUT OBSERVED PROFILES COMES IN A DIRECTORY OR IS A FILE
		if(configCrontrolFile.t1 == 0 && configCrontrolFile.t2 ==0){ // then process only one file
			numberOfFileSpectra = 1;
        	vInputFileSpectra = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
        	strcpy(vInputFileSpectra[0].name,configCrontrolFile.ObservedProfiles);

        	vOutputNameModels = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
        	strcpy(vOutputNameModels[0].name,get_basefilename(configCrontrolFile.InitialGuessModel));	
			//strcat(vOutputNameModels[0].name,"_mod");
			strcat(vOutputNameModels[0].name,MOD_FITS);

			vOutputNameSynthesisAdjusted = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
			strcpy(vOutputNameSynthesisAdjusted[0].name,get_basefilename(configCrontrolFile.ObservedProfiles));
			strcat(vOutputNameSynthesisAdjusted[0].name,STOKES_FIT_EXT);
			
		}
		else
		{
			
			numberOfFileSpectra = (configCrontrolFile.t2 - configCrontrolFile.t1)+1;
			vInputFileSpectra = (nameFile *) malloc(numberOfFileSpectra*sizeof(nameFile));
			vOutputNameModels = (nameFile *) malloc(numberOfFileSpectra*sizeof(nameFile));
			vOutputNameSynthesisAdjusted = (nameFile *) malloc(numberOfFileSpectra*sizeof(nameFile));
			
			int indexName = 0;

			for(i=configCrontrolFile.t1;i<=configCrontrolFile.t2;i++){
				char strIndex[5];
				if(i>=0 && i<10)
					sprintf(strIndex, "0%d", i);
				else
					sprintf(strIndex, "%d", i);
				// FILE NAMES FOR INPUT IMAGES
				strcpy(vInputFileSpectra[indexName].name, configCrontrolFile.ObservedProfiles);
				strcat(vInputFileSpectra[indexName].name,strIndex);
				strcat(vInputFileSpectra[indexName].name,FITS_FILE);
				// FILE NAME FOR OUTPUT MODELS 
				strcpy(vOutputNameModels[indexName].name, configCrontrolFile.ObservedProfiles);
				strcat(vOutputNameModels[indexName].name, strIndex);
				strcat(vOutputNameModels[indexName].name, "_mod");
				if(configCrontrolFile.outputPrefix[0]!='\0'){
					strcat(vOutputNameModels[indexName].name, "_");
					strcat(vOutputNameModels[indexName].name, configCrontrolFile.outputPrefix);
				}
				strcat(vOutputNameModels[indexName].name,FITS_FILE);
				// FILE NAME FOR ADJUSTED SYNTHESIS 
				strcpy(vOutputNameSynthesisAdjusted[indexName].name, configCrontrolFile.ObservedProfiles);
				strcat(vOutputNameSynthesisAdjusted[indexName].name, strIndex);
				strcat(vOutputNameSynthesisAdjusted[indexName].name, "_stokes");
				if(configCrontrolFile.outputPrefix[0]!='\0'){
					strcat(vOutputNameSynthesisAdjusted[indexName].name, "_");
					strcat(vOutputNameSynthesisAdjusted[indexName].name, configCrontrolFile.outputPrefix);
				}
				strcat(vOutputNameSynthesisAdjusted[indexName].name,FITS_FILE);

				indexName++;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
	//  BROADCAST THE NUMBER OF LAMBDAS READS FROM THE FILE AND THE NUMBER OF PIXELS
	MPI_Bcast(&numberOfFileSpectra, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	/**
	 * What we do know is the following: 
	 * --- If the number of files is greater than number of available processors , then scatter the list of files proporcionaly between N processores. 
	 * --- If the number of files is less than number of available processors, then we process each fits file sequentially doing a scatter of N pixels in the each Fits File. 
	 * */

	
	int numFilesPerProcess = numberOfFileSpectra / numProcs;
	int numFilesPerProcessParallel = numberOfFileSpectra % numProcs;
	int sum = 0;                // Sum of counts. Used to calculate displacements

	for ( i = 0; i < numProcs; i++) {
		sendcountsNameInputFiles[i] = numFilesPerProcess;
		displsNameInputFiles[i] = sum;
		sum += sendcountsNameInputFiles[i];
	}

	//printf("\n IDPROC: %d NÚMERO DE FICHEROS POR PROCESO PARALELO %d , NUMERO DE FICHEROS POR PROCESO %d\n",idProc,numFilesPerProcessParallel,numFilesPerProcess);

	
	if(idProc == root && numFilesPerProcess>=1){
		vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
		vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
		vOutputNameSynthesisAdjustedParallel = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
		
		for(i=0;i<numFilesPerProcessParallel;i++){
			strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[(numFilesPerProcess*numProcs)+i].name);
			strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[(numFilesPerProcess*numProcs)+i].name);
			strcpy(vOutputNameSynthesisAdjustedParallel[i].name,vOutputNameSynthesisAdjusted[(numFilesPerProcess*numProcs)+i].name);
		}
		nameFile * auxInput  = (nameFile *)malloc((numFilesPerProcess*numProcs)*sizeof(nameFile));
		nameFile * auxOutput = (nameFile *)malloc((numFilesPerProcess*numProcs)*sizeof(nameFile));
		nameFile * auxOutputSynthesisAdjusted = (nameFile *)malloc((numFilesPerProcess*numProcs)*sizeof(nameFile));

		for(i=0;i<(numFilesPerProcess*numProcs);i++){
			strcpy(auxInput[i].name,vInputFileSpectra[i].name);
			strcpy(auxOutput[i].name,vOutputNameModels[i].name);
			strcpy(auxOutputSynthesisAdjusted[i].name,vOutputNameSynthesisAdjusted[i].name);
		}		
		free(vInputFileSpectra);
		free(vOutputNameModels);
		free(vOutputNameSynthesisAdjusted);
		vInputFileSpectra = auxInput;
		vOutputNameModels = auxOutput;
		vOutputNameSynthesisAdjusted = auxOutputSynthesisAdjusted;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(numberOfFileSpectra>numProcs){ // Scatter name of files 
		Init_Model *vModels;

		vInputFileSpectraLocal = (nameFile *) malloc(sendcountsNameInputFiles[idProc]*sizeof(nameFile));
		vOutputNameModelsLocal = (nameFile *) malloc(sendcountsNameInputFiles[idProc]*sizeof(nameFile));
		vOutputNameSynthesisAdjustedLocal = (nameFile *) malloc(sendcountsNameInputFiles[idProc]*sizeof(nameFile));
		
		if( root == idProc){
			MPI_Scatterv(vInputFileSpectra, sendcountsNameInputFiles, displsNameInputFiles, mpiName, vInputFileSpectraLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(vOutputNameModels, sendcountsNameInputFiles, displsNameInputFiles, mpiName, vOutputNameModelsLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(vOutputNameSynthesisAdjusted, sendcountsNameInputFiles, displsNameInputFiles, mpiName, vOutputNameSynthesisAdjustedLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
		}
		else{
			MPI_Scatterv(NULL, NULL,NULL, mpiName, vInputFileSpectraLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(NULL, NULL,NULL, mpiName, vOutputNameModelsLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(NULL, NULL,NULL, mpiName, vOutputNameSynthesisAdjustedLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//  PROCESS INVERSION OVER EACH FILE ON THE CURRENT PROCESSOR 
		int indexInputFits;
		for(indexInputFits=0;indexInputFits<sendcountsNameInputFiles[idProc];indexInputFits++){
			/****************************************************************************************************/
			// READ PIXELS FROM IMAGE 
			PRECISION timeReadImage,timeExecuteClassicalEstimates;
			clock_t t;
			t = clock();
			fitsImage = readFitsSpectroImage(vInputFileSpectraLocal[indexInputFits].name,0);
			t = clock() - t;
			timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			
			printf("\n************************************************************************************************************************");
			printf("\n\n IDPROC: %d -->  TIME TO READ FITS IMAGE (%s):  %f seconds to execute. NUMBERS OF PIXELS READ %d \n",idProc, vInputFileSpectraLocal[indexInputFits].name, timeReadImage,fitsImage->numPixels); 
			printf("\n************************************************************************************************************************");

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
					//imageStokesAdjust->vLambdaImagen = calloc(imageStokesAdjust->numPixels*imageStokesAdjust->nLambdas, sizeof(PRECISION));
					//imageStokesAdjust->spectroImagen = calloc(imageStokesAdjust->numPixels*imageStokesAdjust->nLambdas*imageStokesAdjust->numStokes, sizeof(PRECISION));
					for( i=0;i<imageStokesAdjust->numPixels;i++){
						imageStokesAdjust->pixels[i].spectro = calloc (nlambda*NPARMS,sizeof(float));
						//imageStokesAdjust->pixels[i].vLambda = calloc (nlambda, sizeof(PRECISION));
						imageStokesAdjust->pixels[i].nLambda = nlambda;
					}
				}


				// INTERPOLATE PSF WITH ARRAY OF LAMBDA READ
				/****************************************************************************************************/	
				// THE NUMBER OF LAMBDAS IS READ FROM INPUT FILES 
				//nlambda = fitsImage->nLambdas;
				// check if read stray light
				// COPY LAMBDA READ IN THE TOP OF FILE 
				int contLambda = 0;
				/*for( i=0;i<fitsImage->numPixels;i++){
					for( j=0;j<nlambda;j++){
						fitsImage->pixels[i].vLambda[j]=vGlobalLambda[j];
						fitsImage->vLambdaImagen[contLambda++] = vGlobalLambda[j];
					}
				}*/

				//***************************************** INIT MEMORY WITH SIZE OF LAMBDA ****************************************************//
				InitializePointerShareCalculation();
				AllocateMemoryDerivedSynthesis(nlambda);
				
				int indexPixel = 0;

				// ALLOCATE MEMORY FOR STORE THE RESULTS 

				vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
				vChisqrf = calloc (fitsImage->numPixels , sizeof(PRECISION));
				vNumIter = calloc (fitsImage->numPixels, sizeof(int));
				t = clock();
				
				printf("\n************************************************************************************************************************\n");
				printf("\nIDPROC: %d -->  DOING INVERSION: %s  \n",idProc,vInputFileSpectraLocal[indexInputFits].name);
				printf("\n************************************************************************************************************************");
				//slog_info(0,"\n***********************  DOING INVERSION *******************************\n\n");

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
					
					PRECISION * slightPixel;
					if(slight==NULL) 
						slightPixel = NULL;
					else{
						if(dimStrayLight==nlambda) 
							slightPixel = slight;
						else 
							slightPixel = slight+nlambda*indexPixel;
					}
					/*lm_mils(cuantic, wlines, fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, fitsImage->pixels[indexPixel].nLambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
							configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);						*/
					lm_mils(cuantic, wlines, vGlobalLambda, nlambda, fitsImage->pixels[indexPixel].spectro, nlambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
							configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);						

					vModels[indexPixel] = initModel;
					if(configCrontrolFile.SaveSynthesisAdjusted){
						int kk;
						for (kk = 0; kk < (nlambda * NPARMS); kk++)
						{
							imageStokesAdjust->pixels[indexPixel].spectro[kk] = spectra[kk] ;
						}						
					}
					//vChisqrf[indexPixel] = chisqrf;
					
					//printf ("\t\t %.2f seconds -- %.2f %%\r",  ((PRECISION)(clock() - t)/CLOCKS_PER_SEC) , ((indexPixel*100.)/fitsImage->numPixels));
				}
				t = clock() - t;
				
				timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
				
				if(!writeFitsImageModels(vOutputNameModelsLocal[indexInputFits].name,fitsImage->rows,fitsImage->cols,vModels,vChisqrf,vNumIter,configCrontrolFile.saveChisqr)){
					printf("\n ERROR WRITING FILE OF MODELS: %s",vOutputNameModelsLocal[indexInputFits].name);
				}

				// PROCESS FILE OF SYNTETIC PROFILES

				if(configCrontrolFile.SaveSynthesisAdjusted){
					// WRITE SINTHETIC PROFILES TO FITS FILE
					if(!writeFitsImageProfiles(vOutputNameSynthesisAdjustedLocal[indexInputFits].name,vInputFileSpectraLocal[indexInputFits].name,imageStokesAdjust)){
						printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",vOutputNameSynthesisAdjustedLocal[indexInputFits].name);
					}
				}
				if(imageStokesAdjust!=NULL){
					for( i=0;i<imageStokesAdjust->numPixels;i++){
						free(imageStokesAdjust->pixels[i].spectro);
					}
					free(imageStokesAdjust->pixels);
					free(imageStokesAdjust);
				}
				free(vModels);
				free(vChisqrf);
				free(vNumIter);

				printf("\n************************************************************************************************************************\n");
				printf(" \n IDPROC: %d --> IMAGE INVERSION FOR IMAGE %s ¡¡¡DONE!!!. TIME: %f *********************\n", idProc, vInputFileSpectraLocal[indexInputFits].name,timeReadImage);
				printf("\n************************************************************************************************************************\n");
			}
			else{
				printf("\n\n IDPROC: %d --> FITS FILE: %s WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n",idProc, vInputFileSpectraLocal[indexInputFits].name);
				//slog_error(0,"\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
			}

			



			freeFitsImage(fitsImage);
			FreeMemoryDerivedSynthesis();
		}
		
		free(vInputFileSpectraLocal);
		free(vOutputNameModelsLocal);
		free(vOutputNameSynthesisAdjustedLocal);
	}
	else{ // all files as parallel 
		if(idProc==root){
			if(vInputFileSpectraParalell!=NULL)
				free(vInputFileSpectraParalell);
			if(vOutputNameModelsParalell!=NULL)
				free(vOutputNameModelsParalell);
			if(vOutputNameSynthesisAdjustedParallel!=NULL)
				free(vOutputNameSynthesisAdjustedParallel);
			
			numFilesPerProcessParallel = numberOfFileSpectra;
			vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameSynthesisAdjustedParallel = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));

			for(i=0;i<numFilesPerProcessParallel;i++){
				strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[i].name);
				strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[i].name);
				strcpy(vOutputNameSynthesisAdjustedParallel[i].name,vOutputNameSynthesisAdjusted[i].name);
			}
			free(vInputFileSpectra);
			free(vOutputNameModels);
			free(vOutputNameSynthesisAdjusted);
		}
		MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
		MPI_Bcast(&numFilesPerProcessParallel, 1, MPI_INT, root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);		
	}
	
	

	MPI_Barrier(MPI_COMM_WORLD);
	int indexInputFits;
	for(indexInputFits=0;indexInputFits<numFilesPerProcessParallel;indexInputFits++){
			//  IF NOT EMPTY READ stray light file
		numPixels=0;
		if(idProc==root){

			if((access(vInputFileSpectraParalell[indexInputFits].name,F_OK)!=-1)){
				
				clock_t t = clock();
				fitsImage = readFitsSpectroImage(vInputFileSpectraParalell[indexInputFits].name,1);
				t = clock() - t;
				PRECISION timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n TIME TO READ FITS IMAGE %s:  %f seconds to execute . NUMBER OF PIXELS READ: %d \n",vInputFileSpectraParalell[indexInputFits].name, timeReadImage,fitsImage->numPixels); 
				
				// COPY LAMBDA READ IN THE TOP OF FILE 
				int contLambda = 0;
				/*for( i=0;i<fitsImage->numPixels;i++){
					for( j=0;j<nlambda;j++){
						fitsImage->pixels[i].vLambda[j]=vGlobalLambda[j];
						fitsImage->vLambdaImagen[contLambda++] = vGlobalLambda[j];
					}
				}*/
				numPixels = fitsImage->numPixels;
			}
			
		}

		MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
		//  BROADCAST THE NUMBER OF PIXELS
		MPI_Bcast(&numPixels, 1, MPI_INT, root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD); // WAIT UNTIL G HAS BEEN READ

		// IF THE NUMBER OF PIXELS IS NOT GREATER THAN 0 WE DON'T CONITUNUE 

		if(numPixels > 0){

			if(idProc == root){
				printf("\n***********************  DOING INVERSION: %s *******************************\n\n",vInputFileSpectraParalell[indexInputFits].name );
				resultsInitModelTotal = calloc (numPixels , sizeof(Init_Model));
				chisqrfTotal = calloc (numPixels , sizeof(PRECISION));
				vNumIterTotal = calloc (numPixels, sizeof(int));
				vSpectraAjustedTotal = calloc (numPixels*nlambda*NPARMS,sizeof(float));
			}
			// allocate memory in all processes 
			InitializePointerShareCalculation();
			AllocateMemoryDerivedSynthesis(nlambda);
			
			
			int numPixelsProceso = numPixels/numProcs;
			int resto = numPixels % numProcs;
			int sum = 0;                // Sum of counts. Used to calculate displacements
			int sumSpectro = 0;
			int sumLambda = 0;
			

			for ( i = 0; i < numProcs; i++) {
				sendcountsPixels[i] = numPixelsProceso;
				if (resto > 0) {
						sendcountsPixels[i]++;
						resto--;
				}
				sendcountsSpectro[i] = sendcountsPixels[i]*nlambda*NPARMS;
				sendcountsLambda[i] = sendcountsPixels[i]*nlambda;
				displsPixels[i] = sum;
				displsSpectro[i] = sumSpectro;
				displsLambda[i] = sumLambda;
				sum += sendcountsPixels[i];
				sumSpectro += sendcountsSpectro[i];
				sumLambda += sendcountsLambda[i];
			}

			MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their vlambda
			local_start = MPI_Wtime();

			// SCATTER VPIXELS 
			vSpectraSplit = calloc(sendcountsSpectro[idProc],sizeof(float));
			//vLambdaSplit = calloc(sendcountsLambda[idProc],sizeof(PRECISION));
			if(configCrontrolFile.SaveSynthesisAdjusted)
				vSpectraAdjustedSplit = calloc(sendcountsSpectro[idProc],sizeof(float));
			
			local_start_scatter = MPI_Wtime();
			
			if( root == idProc){
				MPI_Scatterv(fitsImage->spectroImagen, sendcountsSpectro, displsSpectro, MPI_FLOAT, vSpectraSplit, sendcountsSpectro[idProc], MPI_FLOAT, root, MPI_COMM_WORLD);
				//MPI_Scatterv(fitsImage->vLambdaImagen, sendcountsLambda, displsLambda, MPI_DOUBLE,vLambdaSplit, sendcountsLambda[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD);
			}
			else{
				MPI_Scatterv(NULL, NULL,NULL, MPI_DOUBLE, vSpectraSplit, sendcountsSpectro[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD);
				//MPI_Scatterv(NULL, NULL,NULL, MPI_DOUBLE, vLambdaSplit, sendcountsLambda[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD);
			}		
			local_finish_scatter = MPI_Wtime();

			resultsInitModel = calloc(sendcountsPixels[idProc], sizeof(Init_Model));
			vChisqrf = calloc(sendcountsPixels[idProc], sizeof(PRECISION));
			vNumIter = calloc(sendcountsPixels[idProc], sizeof(int));
			

			local_start_execution = MPI_Wtime();
			for(indexPixel = 0; indexPixel < sendcountsPixels[idProc]; indexPixel++){
				//Initial Model
				Init_Model initModel;
				initModel.eta0 = INITIAL_MODEL.eta0;
				initModel.B = INITIAL_MODEL.B; 
				initModel.gm = INITIAL_MODEL.gm;
				initModel.az = INITIAL_MODEL.az;
				initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
				initModel.mac = INITIAL_MODEL.mac;
				initModel.dopp = INITIAL_MODEL.dopp;
				initModel.aa = INITIAL_MODEL.aa;
				initModel.alfa = INITIAL_MODEL.alfa; 
				initModel.S0 = INITIAL_MODEL.S0;
				initModel.S1 = INITIAL_MODEL.S1;

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
				/*lm_mils(cuantic, wlines, vLambdaSplit+(indexPixel*(nlambda)), nlambda, vSpectraSplit+(indexPixel*(nlambda*NPARMS)), nlambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
					configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);																		*/
				lm_mils(cuantic, wlines, vGlobalLambda, nlambda, vSpectraSplit+(indexPixel*(nlambda*NPARMS)), nlambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
					configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);																							
				
				resultsInitModel[indexPixel] = initModel;
				if(configCrontrolFile.SaveSynthesisAdjusted){
					int kk;
					for (kk = 0; kk < (nlambda * NPARMS); kk++)
					{
						vSpectraAdjustedSplit[ (indexPixel*(nlambda * NPARMS))+kk] = spectra[kk] ;
					}						
				}
				
			}
			local_finish_execution = MPI_Wtime();
			local_start_gather = MPI_Wtime();
			MPI_Igatherv(resultsInitModel, sendcountsPixels[idProc], mpiInitModel, resultsInitModelTotal, sendcountsPixels, displsPixels, mpiInitModel, root, MPI_COMM_WORLD,&mpiRequestSpectro);
			MPI_Igatherv(vChisqrf, sendcountsPixels[idProc], MPI_DOUBLE, chisqrfTotal, sendcountsPixels, displsPixels, MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestLambda);		
			MPI_Gatherv(vNumIter, sendcountsPixels[idProc], MPI_INT, vNumIterTotal, sendcountsPixels, displsPixels, MPI_INT, root, MPI_COMM_WORLD);		
			
			if(configCrontrolFile.SaveSynthesisAdjusted)
				MPI_Gatherv(vSpectraAdjustedSplit, sendcountsSpectro[idProc], MPI_FLOAT, vSpectraAjustedTotal, sendcountsSpectro, displsSpectro, MPI_FLOAT, root, MPI_COMM_WORLD);		
				
			MPI_Wait(&mpiRequestSpectro, MPI_STATUS_IGNORE);
			MPI_Wait(&mpiRequestLambda, MPI_STATUS_IGNORE);		
			
			local_finish_gather = MPI_Wtime();
			local_finish = MPI_Wtime();
			local_elapsed = local_finish - local_start;
			local_elapsed_execution = local_finish_execution - local_start_execution;
			local_elapsed_scatter = local_finish_scatter - local_start_scatter;
			local_elapsed_gather = local_finish_gather - local_start_gather;
			MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&local_elapsed_execution, &elapsed_execution, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&local_elapsed_scatter, &elapsed_scatter, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&local_elapsed_gather, &elapsed_gather, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			if(idProc==root){		 
				printf("\n Elapsed SCATTER time = %lf seconds\n", elapsed_scatter);
				printf("\n***********************************\n");
				printf("\n Elapsed GATHER time = %lf seconds\n", elapsed_gather);
				printf("\n***********************************\n");			
				printf("\n Elapsed EXECUTION time = %lf seconds\n", elapsed_execution);
				printf("\n***********************************\n");
				printf("\n Elapsed TOTAL time = %lf seconds\n", elapsed);
				printf("\n***********************************\n");
				double timeWriteImage;
				clock_t t;
				t = clock();

				if(!writeFitsImageModels(vOutputNameModelsParalell[indexInputFits].name,fitsImage->rows,fitsImage->cols,resultsInitModelTotal,chisqrfTotal,vNumIterTotal,configCrontrolFile.saveChisqr)){
						printf("\n ERROR WRITING FILE OF MODELS: %s",vOutputNameModelsParalell[indexInputFits].name);
				}
				t = clock() - t;
				timeWriteImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n TIME TO WRITE FITS IMAGE:  %f seconds to execute \n", timeWriteImage); 
				// PROCESS FILE OF SYNTETIC PROFILES
				if(configCrontrolFile.SaveSynthesisAdjusted){

					for(indexPixel=0;indexPixel<numPixels;indexPixel++)
					{	
						int kk;
						for (kk = 0; kk < (nlambda * NPARMS); kk++)
						{
							fitsImage->pixels[indexPixel].spectro[kk] = vSpectraAjustedTotal[kk+(indexPixel*(nlambda * NPARMS))] ;
						}
					}					
					// WRITE SINTHETIC PROFILES TO FITS FILE
					if(!writeFitsImageProfiles(vOutputNameSynthesisAdjustedParallel[indexInputFits].name,vInputFileSpectraParalell[indexInputFits].name,fitsImage)){
						printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",vOutputNameSynthesisAdjustedParallel[indexInputFits].name);
					}
				}

				free(resultsInitModelTotal);		
				free(chisqrfTotal);
				free(vNumIterTotal);
				if(configCrontrolFile.SaveSynthesisAdjusted){
					free(vSpectraAjustedTotal);
				}
				printf("\n************************************************************************************************************************\n");
				printf("\n INVERSION OF IMAGE %s ¡¡¡DONE!!!\n", vInputFileSpectraParalell[indexInputFits].name);
				printf("\n************************************************************************************************************************\n");
			}
			else{
				if(configCrontrolFile.SaveSynthesisAdjusted)
					free(vSpectraAdjustedSplit);
				free(vSpectraSplit);
				//free(vLambdaSplit);
				free(resultsInitModel);				
				free(vChisqrf);
				free(vNumIter);
			}
			
			FreeMemoryDerivedSynthesis();
			
		}
		else{
			if(idProc==root){
				printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT %s ******************************",vInputFileSpectraParalell[indexInputFits].name);
			}
		}
		if(idProc==root){
			freeFitsImage(fitsImage);
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
	//printf("\n\n TOTAL sec : %.16g segundos\n", total_secs);
	free(cuantic);
	free(wlines);
	free(vGlobalLambda);
	//FreeMemoryDerivedSynthesis();
	// FREE TYPE OF MPI
	MPI_Type_free(&mpiInitModel);
	MPI_Finalize() ;
	free(G);
	return 0;
}

