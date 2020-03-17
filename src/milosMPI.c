
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
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


// ***************************** FUNCTIONS TO READ FITS FILE *********************************************************


Cuantic *cuantic; // Global variable with cuantic information 


REAL *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
REAL *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
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


PRECISION *GMAC,*GMAC_DERIV, *G; // GAUSSIAN MUST BE IN DOUBLE PRECISION 
PRECISION *dirConvPar; // AUX GLOBAL VECTOR for calculate direct convolutions
//REAL *resultConv; // aux global vector for store direct convolution
REAL AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
REAL * opa;
int FGlobal, HGlobal, uuGlobal;

//PRECISION *d_spectra, *spectra, *spectra_mac;
REAL *d_spectra, *spectra, *spectra_mac;


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

_Complex double  *z,* zden, * zdiv;

gsl_vector *eval;
gsl_matrix *evec;
gsl_eigen_symmv_workspace * workspace;

int main(int argc, char **argv)
{
	int i;  // indexes 
	
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

	/************************************/
	const int root=0;	
	
	// FINISH STARTING PROGRAM 

	PRECISION *wlines;
	int nlambda, numPixels, indexPixel;
	
	// allocate memory for eigen values
	eval = gsl_vector_alloc (NTERMS);
  	evec = gsl_matrix_alloc (NTERMS, NTERMS);
	workspace = gsl_eigen_symmv_alloc (NTERMS);

	//*****
	Init_Model INITIAL_MODEL;
	PRECISION * deltaLambda, * PSF;
	int N_SAMPLES_PSF;
	
	
	// CONFIGURACION DE PARAMETROS A INVERTIR
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	//----------------------------------------------

	PRECISION * slight = NULL;
	
	int dimStrayLight;
	int numberOfFileSpectra;

   	nameFile * vInputFileSpectra = NULL;
	nameFile * vInputFileSpectraParalell = NULL;
	nameFile * vInputFileSpectraDiv2Parallel = NULL;
	nameFile * vOutputNameModels = NULL;
	nameFile * vOutputNameModelsParalell = NULL;
	nameFile * vOutputNameModelsDiv2Parallel = NULL;
	nameFile * vOutputNameSynthesisAdjusted = NULL;
	nameFile * vOutputNameSynthesisAdjustedParallel = NULL;
	nameFile * vOutputNameSynthesisAdjustedDiv2Parallel = NULL;
	nameFile * vInputFileSpectraLocal = NULL;
	nameFile * vOutputNameModelsLocal = NULL;
	nameFile * vOutputNameSynthesisAdjustedLocal = NULL;

	
	const char	* nameInputFilePSF ;

   FitsImage * fitsImage = NULL;
	PRECISION  dat[7];

	double local_start, local_finish, local_elapsed, elapsed;
	double local_start_execution, local_finish_execution, local_elapsed_execution, elapsed_execution;
	double local_start_scatter, local_finish_scatter, local_elapsed_scatter, elapsed_scatter;
	double local_start_gather, local_finish_gather, local_elapsed_gather, elapsed_gather;
	
	Init_Model * resultsInitModel;
	Init_Model * resultsInitModelTotal;
	float * chisqrfTotal, *vChisqrf;
	int * vNumIter, * vNumIterTotal; // to store the number of iterations used to converge for each pixel
	float  * vSpectraSplit, * vSpectraAdjustedSplit, * vSpectraAjustedTotal;

	int sendcountsPixels [numProcs] ; // array describing how many elements to send to each process
	int sendcountsSpectro [numProcs];
	int sendcountsLambda [numProcs];
	int sendcountsNameInputFiles [numProcs];  // how many files per process
	int displsPixels [numProcs];  // array describing the displacements where each segment begins
	int displsSpectro [numProcs];
	//int displsLambda [numProcs];
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
	PRECISION * vGlobalLambda, *vOffsetsLambda;

	if(configCrontrolFile.useMallaGrid){ // read lambda from grid file
		indexLine = readMallaGrid(configCrontrolFile.MallaGrid, &initialLambda, &step, &finalLambda, (idProc==root));      
		nlambda = ((finalLambda-initialLambda)/step)+1;
		vOffsetsLambda = calloc(nlambda,sizeof(PRECISION));
		vOffsetsLambda[0] = initialLambda;
		for(i=1;i<nlambda;i++){
			vOffsetsLambda[i] = vOffsetsLambda[i-1]+step;
		}
		// pass to armstrong 
		initialLambda = initialLambda/1000;
		step = step/1000;
		finalLambda = finalLambda/1000;
		vGlobalLambda = calloc(nlambda,sizeof(PRECISION));
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(configCrontrolFile.AtomicParametersFile,dat,indexLine,(idProc==root));
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n QUANTUM LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVELENGHT: %f",configCrontrolFile.CentralWaveLenght);
			exit(1);
		}
		vGlobalLambda[0]=configCrontrolFile.CentralWaveLenght+(initialLambda);
   		for(i=1;i<nlambda;i++){
        	vGlobalLambda[i]=vGlobalLambda[i-1]+step;
     	}
	}
	else{
		vGlobalLambda = readFitsLambdaToArray(configCrontrolFile.WavelengthFile,&indexLine,&nlambda);
		if(vGlobalLambda==NULL){
			printf("\n FILE WITH WAVELENGHT HAS NOT BEEN READ PROPERLY, please check it.\n");
			free(vGlobalLambda);
			exit(EXIT_FAILURE);
		}
		configCrontrolFile.CentralWaveLenght = readFileCuanticLines(configCrontrolFile.AtomicParametersFile,dat,indexLine,(idProc==root));
		if(configCrontrolFile.CentralWaveLenght==0){
			printf("\n QUANTUM LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
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
	

	/**************************************** READ FITS  STRAY LIGHT ******************************/
	MPI_Barrier(MPI_COMM_WORLD);
	if(access(configCrontrolFile.StrayLightFile,F_OK)!=-1){ //  IF NOT EMPTY READ stray light file 
		slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,nlambda);
	}

	
	// ************************** DEFINE PLANS TO EXECUTE MACROTURBULENCE IF NECESSARY **********************************************//
	// MACROTURBULENCE PLANS
	
	inFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC = fftw_plan_dft_1d(nlambda, inFilterMAC, outFilterMAC, FFT_FORWARD, FFTW_MEASURE     );
	inFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC_DERIV = fftw_plan_dft_1d(nlambda, inFilterMAC_DERIV, outFilterMAC_DERIV, FFT_FORWARD, FFTW_MEASURE     );


	inSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planForwardMAC = fftw_plan_dft_1d(nlambda, inSpectraFwMAC, outSpectraFwMAC, FFT_FORWARD, FFTW_MEASURE    );
	inSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
	planBackwardMAC = fftw_plan_dft_1d(nlambda, inSpectraBwMAC, outSpectraBwMAC, FFT_BACKWARD, FFTW_MEASURE    );

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
					readPSFFile(deltaLambda,PSF,nameInputFilePSF,configCrontrolFile.CentralWaveLenght);
					// CHECK if values of deltaLambda are in the same range of vLambda. For do that we truncate to 4 decimal places 
					if( (trunc(vOffsetsLambda[0])) < (trunc(deltaLambda[0]))  || (trunc(vOffsetsLambda[nlambda-1])) > (trunc(deltaLambda[N_SAMPLES_PSF-1])) ){
						if(idProc==root){
							printf("\n\n ERROR: The wavelength range given in the PSF file is smaller than the range in the mesh file [%lf,%lf] [%lf,%lf]  \n\n",deltaLambda[0],vOffsetsLambda[0],deltaLambda[N_SAMPLES_PSF-1],vOffsetsLambda[nlambda-1]);
						}
						exit(EXIT_FAILURE);
					}
					G = calloc(nlambda,sizeof(PRECISION));
					double offset =0;
					int posWL = 0;
					for(i=0;i<nlambda && !posWL;i++){
						//if( (trunc(vGlobalLambda[i]*1000)/1000)== (trunc(configCrontrolFile.CentralWaveLenght*1000)/1000))
						if( fabs(trunc(vOffsetsLambda[i]))==0)
							posWL = i;
					}
					if(posWL!= (nlambda/2)){ // move center to the middle of samples
						//printf("\nPOS CENTRAL WL %i",posWL);
						offset = (((nlambda/2)-posWL)*step)*1000;
						//printf ("\n OFFSET IS %f\n",offset);
					}					
					interpolationLinearPSF(deltaLambda,  PSF, vOffsetsLambda , N_SAMPLES_PSF, G, nlambda,offset);		
					sizeG=nlambda;	
				}
				else{
					//G = fgauss_WL(FWHM,vGlobalLambda[1]-vGlobalLambda[0],vGlobalLambda[0],vGlobalLambda[nlambda/2],nlambda,&sizeG);
					if(idProc==root)
						printf("\n****************** ERROR THE PSF FILE is empty or damaged.******************\n");
					exit(EXIT_FAILURE);
				}
		}

		
		//PSF FILTER PLANS 
		inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF = fftw_plan_dft_1d(nlambda, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_MEASURE     );
		inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
		planBackwardPSF = fftw_plan_dft_1d(nlambda, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_MEASURE     );

		fftw_complex * in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		int i;
		for (i = 0; i < nlambda; i++)
		{
			in[i] = G[i] + 0 * _Complex_I;
		}
		fftw_G_PSF = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_plan p = fftw_plan_dft_1d(nlambda, in, fftw_G_PSF, FFT_FORWARD, FFTW_MEASURE   );
		fftw_execute(p);
		for (i = 0; i < nlambda; i++)
		{
			fftw_G_PSF[i] = fftw_G_PSF[i] / nlambda;
		}
		fftw_destroy_plan(p);
		fftw_free(in);
		
		inPSF_MAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC = fftw_plan_dft_1d(nlambda, inPSF_MAC, fftw_G_MAC_PSF, FFT_FORWARD, FFTW_MEASURE     );
		inMulMacPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFilters = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC = fftw_plan_dft_1d(nlambda, inMulMacPSF, outConvFilters, FFT_BACKWARD, FFTW_MEASURE     );


		inPSF_MAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_DERIV_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inPSF_MAC_DERIV, fftw_G_MAC_DERIV_PSF, FFT_FORWARD, FFTW_MEASURE     );
		inMulMacPSFDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFiltersDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inMulMacPSFDeriv, outConvFiltersDeriv, FFT_BACKWARD, FFTW_MEASURE     );			

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
	int numFilesPer2ProcessParallel=0;
	if(numFilesPerProcessParallel>=(numProcs/2)){ // DIVIDE EACH FILE IN TWO PROCESS
		numFilesPer2ProcessParallel = numProcs/2;
		numFilesPerProcessParallel = numFilesPerProcessParallel - (numProcs/2);
	}
	else
	{
		numFilesPer2ProcessParallel = 0;
	}

	int sum = 0;                // Sum of counts. Used to calculate displacements

	for ( i = 0; i < numProcs; i++) {
		sendcountsNameInputFiles[i] = numFilesPerProcess;
		displsNameInputFiles[i] = sum;
		sum += sendcountsNameInputFiles[i];
	}

	//**************************************** CREATE GROUPS FOR DIVIDE IMAGE IN 2 ********************************************/
	
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	int numGroups = numProcs/2;
	MPI_Group vGroups [numGroups]; // will contain the numbers of process of each group

	if(numProcs%2 ==0){
		for(i=0;i<numProcs;i=i+2){
			int ranks[2];
			ranks[0] = i;
			ranks[1] = i+1;
			MPI_Group_incl(world_group, 2, ranks, &vGroups[i/2]);
		}
	}
	else{
		int indProc=0;
		for(i=0;i<numGroups;i++){
			if(i==(numGroups-1)){ // add resto
				int ranks[3];
				ranks[0]=indProc++;
				ranks[1]=indProc++;
				ranks[2]=indProc++;
				MPI_Group_incl(world_group, 3, ranks, &vGroups[i]);
			}
			else{
				int ranks[2];
				ranks[0]=indProc++;
				ranks[1]=indProc++;
				MPI_Group_incl(world_group, 2, ranks, &vGroups[i]);
			}
		}
	}
	// Create the new communicator from that group of processes.
	MPI_Comm vCommunicators[numGroups];
	for(i=0;i<numGroups;i++){
		MPI_Comm_create(MPI_COMM_WORLD, vGroups[i], &vCommunicators[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	int myGroup = idProc/2;
	if(myGroup==numGroups)
		myGroup = myGroup-1;
	
	int myGroupRank;
	int groupRoot = 0; // process 0 of group will be the root 
	int myGroupSize;
	MPI_Group_rank(vGroups[myGroup], &myGroupRank);		
	MPI_Group_size(vGroups[myGroup], &myGroupSize);
	MPI_Barrier(MPI_COMM_WORLD);
	

	//**************************************** END OF CREATE GROUPS FOR DIVIDE IMAGE IN 2 ********************************************/

	if(idProc == root){
		if(numFilesPerProcess>=1){
			vInputFileSpectraDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
			vOutputNameModelsDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
			vOutputNameSynthesisAdjustedDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
			
			for(i=0;i<numFilesPer2ProcessParallel;i++){
				strcpy(vInputFileSpectraDiv2Parallel[i].name,vInputFileSpectra[(numFilesPerProcess*numProcs)+i].name);
				strcpy(vOutputNameModelsDiv2Parallel[i].name,vOutputNameModels[(numFilesPerProcess*numProcs)+i].name);
				strcpy(vOutputNameSynthesisAdjustedDiv2Parallel[i].name,vOutputNameSynthesisAdjusted[(numFilesPerProcess*numProcs)+i].name);
			}

			if(numFilesPerProcessParallel>0){
				vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
				vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
				vOutputNameSynthesisAdjustedParallel = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
				
				for(i=0;i<numFilesPerProcessParallel;i++){
					strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[(numFilesPerProcess*numProcs)+numFilesPer2ProcessParallel+i].name);
					strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[(numFilesPerProcess*numProcs)+numFilesPer2ProcessParallel+i].name);
					strcpy(vOutputNameSynthesisAdjustedParallel[i].name,vOutputNameSynthesisAdjusted[(numFilesPerProcess*numProcs)+numFilesPer2ProcessParallel+i].name);
				}
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
		else{
			if(vInputFileSpectraDiv2Parallel!=NULL)
				free(vInputFileSpectraDiv2Parallel);
			if(vOutputNameModelsDiv2Parallel!=NULL)
				free(vOutputNameModelsDiv2Parallel);
			if(vOutputNameSynthesisAdjustedDiv2Parallel!=NULL)
				free(vOutputNameSynthesisAdjustedDiv2Parallel);			
			if(vInputFileSpectraParalell!=NULL)
				free(vInputFileSpectraParalell);
			if(vOutputNameModelsParalell!=NULL)
				free(vOutputNameModelsParalell);
			if(vOutputNameSynthesisAdjustedParallel!=NULL)
				free(vOutputNameSynthesisAdjustedParallel);

			if(numFilesPer2ProcessParallel>0){
				vInputFileSpectraDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
				vOutputNameModelsDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
				vOutputNameSynthesisAdjustedDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));				
				for(i=0;i<numFilesPer2ProcessParallel;i++){
					strcpy(vInputFileSpectraDiv2Parallel[i].name,vInputFileSpectra[i].name);
					strcpy(vOutputNameModelsDiv2Parallel[i].name,vOutputNameModels[i].name);
					strcpy(vOutputNameSynthesisAdjustedDiv2Parallel[i].name,vOutputNameSynthesisAdjusted[i].name);
				}
			}
			
			vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameSynthesisAdjustedParallel = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));

			for(i=0;i<numFilesPerProcessParallel;i++){
				strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[numFilesPer2ProcessParallel+i].name);
				strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[numFilesPer2ProcessParallel+i].name);
				strcpy(vOutputNameSynthesisAdjustedParallel[i].name,vOutputNameSynthesisAdjusted[numFilesPer2ProcessParallel+i].name);
			}
			free(vInputFileSpectra);
			vInputFileSpectra = NULL;
			free(vOutputNameModels);
			vOutputNameModels = NULL;
			free(vOutputNameSynthesisAdjusted);
			vOutputNameSynthesisAdjusted = NULL;			
		}
	}

	
	MPI_Barrier(MPI_COMM_WORLD);

	if(numFilesPer2ProcessParallel>0){ // CASE DIVIDE EACH IMAGE BETWEEN TWO PROCESS
		if(idProc!=root){
			vInputFileSpectraDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
			vOutputNameModelsDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));
			vOutputNameSynthesisAdjustedDiv2Parallel = (nameFile *)malloc(numFilesPer2ProcessParallel*sizeof(nameFile));		
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(vInputFileSpectraDiv2Parallel, numFilesPer2ProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Bcast(vOutputNameModelsDiv2Parallel, numFilesPer2ProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Bcast(vOutputNameSynthesisAdjustedDiv2Parallel, numFilesPer2ProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);	
	}

	if(numFilesPerProcessParallel>0){
		if(idProc!=root){
			vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameSynthesisAdjustedParallel = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));		
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(vInputFileSpectraParalell, numFilesPerProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Bcast(vOutputNameModelsParalell, numFilesPerProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Bcast(vOutputNameSynthesisAdjustedParallel, numFilesPerProcessParallel,mpiName , root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);		
	}

	
	// MEMORY IF SCATTER IMAGES 

	FitsImage ** fitsImages;
	fitsImages = (FitsImage **)malloc(numFilesPerProcessParallel*sizeof(FitsImage*));
	MPI_Request * vMpiRequestScatter = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	MPI_Request * vMpiRequestInitModel = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	MPI_Request * vMpiRequestChisqr = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	MPI_Request * vMpiRequestIter = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	MPI_Request * vMpiRequestSpectraAd = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	MPI_Request * vMpiRequestReduceExecution = malloc(numFilesPerProcessParallel*sizeof(MPI_Request));
	int * vNumPixelsImage = malloc(numFilesPerProcessParallel*sizeof(int));
	Init_Model ** resultsInitModel_L = (Init_Model **)malloc(numFilesPerProcessParallel*sizeof(Init_Model*));
	Init_Model ** resultsInitModelTotal_L = (Init_Model **)malloc(numFilesPerProcessParallel*sizeof(Init_Model*));
	float ** chisqrfTotal_L = (float **) malloc(numFilesPerProcessParallel*sizeof(float*));
	float **vChisqrf_L = (float **) malloc(numFilesPerProcessParallel*sizeof(float*));
	int **vNumIter_L = (int **) malloc(numFilesPerProcessParallel*sizeof(int*));
	int **vNumIterTotal_L = (int **) malloc(numFilesPerProcessParallel*sizeof(int*));
	float  **vSpectraSplit_L = (float **) malloc(numFilesPerProcessParallel*sizeof(float*));
	float  **vSpectraAdjustedSplit_L = (float **) malloc(numFilesPerProcessParallel*sizeof(float*));
	float  **vSpectraAjustedTotal_L = (float **) malloc(numFilesPerProcessParallel*sizeof(float*));
	int sendcountsPixels_L [numFilesPerProcessParallel][numProcs] ; // array describing how many elements to send to each process
	int sendcountsSpectro_L [numFilesPerProcessParallel][numProcs];
	int sendcountsLambda_L [numFilesPerProcessParallel][numProcs];
	
	int displsPixels_L [numFilesPerProcessParallel][numProcs];  // array describing the displacements where each segment begins
	int displsSpectro_L [numFilesPerProcessParallel][numProcs];
	//int displsLambda [numProcs];
	
	double * vElapsed_execution = calloc(numFilesPerProcessParallel,sizeof(double));
	int indexInputFits;
	if(numFilesPerProcessParallel){
		// FIRST SCATTER ALL PIXELS BETWEEN ALL PROCESS
		clock_t t = clock();
		for(indexInputFits=0;indexInputFits<numFilesPerProcessParallel;indexInputFits++){
			if(idProc==root){

				if((access(vInputFileSpectraParalell[indexInputFits].name,F_OK)!=-1)){
					clock_t t = clock();
					fitsImages[indexInputFits] = readFitsSpectroImage(vInputFileSpectraParalell[indexInputFits].name,1);
					t = clock() - t;
					PRECISION timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
					printf("\n TIME TO READ FITS IMAGE %s:  %f seconds to execute . NUMBER OF PIXELS READ: %d \n",vInputFileSpectraParalell[indexInputFits].name, timeReadImage,fitsImages[indexInputFits]->numPixels); 
					vNumPixelsImage[indexInputFits] = fitsImages[indexInputFits]->numPixels;
				}
				
			}
			MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
			//  BROADCAST THE NUMBER OF PIXELS
			MPI_Bcast(&vNumPixelsImage[indexInputFits], 1, MPI_INT, root , MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD); // WAIT UNTIL G HAS BEEN READ

			// IF THE NUMBER OF PIXELS IS NOT GREATER THAN 0 WE DON'T CONITUNUE 
			if(vNumPixelsImage[indexInputFits] > 0){
				
				if(idProc == root){
					printf("\n***********************  DOING INVERSION: %s *******************************\n\n",vInputFileSpectraParalell[indexInputFits].name );
					resultsInitModelTotal_L[indexInputFits] = calloc (vNumPixelsImage[indexInputFits] , sizeof(Init_Model));
					chisqrfTotal_L[indexInputFits] = calloc (vNumPixelsImage[indexInputFits] , sizeof(float));
					vNumIterTotal_L[indexInputFits] = calloc (vNumPixelsImage[indexInputFits], sizeof(int));
					if(configCrontrolFile.SaveSynthesisAdjusted)
						vSpectraAjustedTotal_L[indexInputFits] = calloc (vNumPixelsImage[indexInputFits]*nlambda*NPARMS,sizeof(float));
				}
				// allocate memory in all processes 

				int numPixelsProceso = vNumPixelsImage[indexInputFits]/(numProcs);
				int resto = vNumPixelsImage[indexInputFits] % (numProcs);
				int sum = 0;                // Sum of counts. Used to calculate displacements
				int sumSpectro = 0;
				int sumLambda = 0;
				sendcountsPixels_L[indexInputFits][0] = 0;
				sendcountsSpectro_L[indexInputFits][0] = 0;
				sendcountsLambda_L[indexInputFits][0] =0;
				displsPixels_L[indexInputFits][0] =0;
				displsSpectro_L[indexInputFits][0] =0;
				for ( i = 0; i < numProcs; i++) {
					sendcountsPixels_L[indexInputFits][i] = numPixelsProceso;
					if (resto > 0) {
							sendcountsPixels_L[indexInputFits][i]++;
							resto--;
					}
					sendcountsSpectro_L[indexInputFits][i] = (sendcountsPixels_L[indexInputFits][i])*nlambda*NPARMS;
					sendcountsLambda_L[indexInputFits][i] = (sendcountsPixels_L[indexInputFits][i])*nlambda;
					displsPixels_L[indexInputFits][i] = sum;
					displsSpectro_L[indexInputFits][i] = sumSpectro;
					//displsLambda[i] = sumLambda;
					sum += sendcountsPixels_L[indexInputFits][i];
					sumSpectro += sendcountsSpectro_L[indexInputFits][i];
					sumLambda += sendcountsLambda_L[indexInputFits][i];
				}

				MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their vlambda
				local_start = MPI_Wtime();

				// SCATTER VPIXELS 
				vSpectraSplit_L[indexInputFits] = calloc(sendcountsSpectro_L[indexInputFits][idProc],sizeof(float));
				if(configCrontrolFile.SaveSynthesisAdjusted)
					vSpectraAdjustedSplit_L[indexInputFits] = calloc(sendcountsSpectro_L[indexInputFits][idProc],sizeof(float));
				
				local_start_scatter = MPI_Wtime();
				

				MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their vlambda				
				if( root == idProc){
					//MPI_Scatterv(fitsImages[indexInputFits]->spectroImagen, sendcountsSpectro_L[indexInputFits], displsSpectro_L[indexInputFits], MPI_FLOAT, vSpectraSplit_L[indexInputFits], sendcountsSpectro_L[indexInputFits][idProc], MPI_FLOAT, root, MPI_COMM_WORLD);
					MPI_Iscatterv(fitsImages[indexInputFits]->spectroImagen, sendcountsSpectro_L[indexInputFits], displsSpectro_L[indexInputFits], MPI_FLOAT, vSpectraSplit_L[indexInputFits], sendcountsSpectro_L[indexInputFits][idProc], MPI_FLOAT, root, MPI_COMM_WORLD,&vMpiRequestScatter[indexInputFits]);
				}
				else{
					//MPI_Scatterv(NULL, NULL,NULL, MPI_FLOAT, vSpectraSplit_L[indexInputFits], sendcountsSpectro_L[indexInputFits][idProc], MPI_FLOAT, root, MPI_COMM_WORLD);
					MPI_Iscatterv(NULL, NULL,NULL, MPI_FLOAT, vSpectraSplit_L[indexInputFits], sendcountsSpectro_L[indexInputFits][idProc], MPI_FLOAT, root, MPI_COMM_WORLD,&vMpiRequestScatter[indexInputFits]);
				}		
				local_finish_scatter = MPI_Wtime();

				resultsInitModel_L[indexInputFits] = calloc(sendcountsPixels_L[indexInputFits][idProc], sizeof(Init_Model));
				vChisqrf_L[indexInputFits] = calloc(sendcountsPixels_L[indexInputFits][idProc], sizeof(float));
				vNumIter_L[indexInputFits] = calloc(sendcountsPixels_L[indexInputFits][idProc], sizeof(int));
			}
			else if (idProc==root){
				printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT %s ******************************",vInputFileSpectraParalell[indexInputFits].name);
			}
		}
		MPI_Waitall(numFilesPerProcessParallel,vMpiRequestScatter,MPI_STATUSES_IGNORE);
		if(idProc==root){
			t = clock() - t;
			PRECISION timeTotalExecution = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			printf("\n TOTAL TIME TO READ AND SCATTER IMAGES: %f",timeTotalExecution);
		}
	}

	AllocateMemoryDerivedSynthesis(nlambda);
	//*************************************** ONE IMAGE PER PROCESSOR *********************************

	if(numFilesPerProcessParallel){

		//clock_t t = clock();
		// EACH PROC PROCESS THER PIXELS 
		for(indexInputFits=0;indexInputFits<numFilesPerProcessParallel;indexInputFits++){
			if(vNumPixelsImage[indexInputFits] > 0){
				local_start_execution = MPI_Wtime();
				for(indexPixel = 0; indexPixel < sendcountsPixels_L[indexInputFits][idProc]; indexPixel++){
					float * vAuxSpectraSplit = vSpectraSplit_L[indexInputFits];
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

					// CLASSICAL ESTIMATES TO GET B, GAMMA, vlos, azimuth
					estimacionesClasicas(wlines[1], vGlobalLambda, nlambda, vAuxSpectraSplit+(indexPixel*(nlambda*NPARMS)), &initModel,1);
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
					
					
					lm_mils(cuantic, wlines, vGlobalLambda, nlambda, vAuxSpectraSplit+(indexPixel*(nlambda*NPARMS)), nlambda, &initModel, spectra, &(vChisqrf_L[indexInputFits][indexPixel]), slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
						configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&(vNumIter_L[indexInputFits][indexPixel]));																		
					
					resultsInitModel_L[indexInputFits][indexPixel] = initModel;
					if(configCrontrolFile.SaveSynthesisAdjusted){
						int kk;
						for (kk = 0; kk < (nlambda * NPARMS); kk++)
						{
							vSpectraAdjustedSplit_L[indexInputFits][ (indexPixel*(nlambda * NPARMS))+kk] = spectra[kk] ;
						}						
					}
					
				}

				MPI_Igatherv(resultsInitModel_L[indexInputFits], sendcountsPixels_L[indexInputFits][idProc], mpiInitModel, resultsInitModelTotal_L[indexInputFits], sendcountsPixels_L[indexInputFits], displsPixels_L[indexInputFits], mpiInitModel, root, MPI_COMM_WORLD,&vMpiRequestInitModel[indexInputFits]);
				MPI_Igatherv(vChisqrf_L[indexInputFits], sendcountsPixels_L[indexInputFits][idProc], MPI_FLOAT, chisqrfTotal_L[indexInputFits], sendcountsPixels_L[indexInputFits], displsPixels_L[indexInputFits], MPI_FLOAT, root, MPI_COMM_WORLD,&vMpiRequestChisqr[indexInputFits]);		
				MPI_Igatherv(vNumIter_L[indexInputFits], sendcountsPixels_L[indexInputFits][idProc], MPI_INT, vNumIterTotal_L[indexInputFits], sendcountsPixels_L[indexInputFits], displsPixels_L[indexInputFits], MPI_INT, root, MPI_COMM_WORLD,&vMpiRequestIter[indexInputFits]);		
				
				if(configCrontrolFile.SaveSynthesisAdjusted)
					MPI_Igatherv(vSpectraAdjustedSplit_L[indexInputFits], sendcountsSpectro_L[indexInputFits][idProc], MPI_FLOAT, vSpectraAjustedTotal_L[indexInputFits], sendcountsSpectro_L[indexInputFits], displsSpectro_L[indexInputFits], MPI_FLOAT, root, MPI_COMM_WORLD,&vMpiRequestSpectraAd[indexInputFits]);		
				local_elapsed_execution = MPI_Wtime() - local_start_execution;
				MPI_Ireduce(&local_elapsed_execution, &vElapsed_execution[indexInputFits], 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD,&vMpiRequestReduceExecution[indexInputFits]);

			}
		}


		if(idProc==root){
			int indexInputFits = 0;
			do{
				MPI_Wait(&vMpiRequestInitModel[indexInputFits],MPI_STATUS_IGNORE);
				MPI_Wait(&vMpiRequestChisqr[indexInputFits],MPI_STATUS_IGNORE);
				MPI_Wait(&vMpiRequestIter[indexInputFits],MPI_STATUS_IGNORE);
				MPI_Wait(&vMpiRequestReduceExecution[indexInputFits],MPI_STATUS_IGNORE);
				if(configCrontrolFile.SaveSynthesisAdjusted)
					MPI_Wait(&vMpiRequestSpectraAd[indexInputFits],MPI_STATUS_IGNORE);

					double timeWriteImage;
					clock_t t;
					t = clock();

					if(!writeFitsImageModels(vOutputNameModelsParalell[indexInputFits].name,fitsImages[indexInputFits]->rows,fitsImages[indexInputFits]->cols,resultsInitModelTotal_L[indexInputFits],chisqrfTotal_L[indexInputFits],vNumIterTotal_L[indexInputFits],configCrontrolFile.saveChisqr)){
							printf("\n ERROR WRITING FILE OF MODELS: %s",vOutputNameModelsParalell[indexInputFits].name);
					}
					t = clock() - t;
					timeWriteImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
					
					// PROCESS FILE OF SYNTETIC PROFILES
					if(configCrontrolFile.SaveSynthesisAdjusted){
						fitsImages[indexInputFits]->pixels = calloc(fitsImages[indexInputFits]->numPixels, sizeof(vpixels));
						for( i=0;i<fitsImages[indexInputFits]->numPixels;i++){
							fitsImages[indexInputFits]->pixels[i].spectro = calloc ((fitsImages[indexInputFits]->numStokes*fitsImages[indexInputFits]->nLambdas),sizeof(float));
							//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
						}		
						for(indexPixel=0;indexPixel<fitsImages[indexInputFits]->numPixels;indexPixel++)
						{	
							int kk;
							for (kk = 0; kk < (nlambda * NPARMS); kk++)
							{
								fitsImages[indexInputFits]->pixels[indexPixel].spectro[kk] = vSpectraAjustedTotal_L[indexInputFits][kk+(indexPixel*(nlambda * NPARMS))] ;
							}
						}					
						// WRITE SINTHETIC PROFILES TO FITS FILE
						if(!writeFitsImageProfiles(vOutputNameSynthesisAdjustedParallel[indexInputFits].name,vInputFileSpectraParalell[indexInputFits].name,fitsImages[indexInputFits])){
							printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",vOutputNameSynthesisAdjustedParallel[indexInputFits].name);
						}
						
						for( i=0;i<fitsImages[indexInputFits]->numPixels;i++){
							free(fitsImages[indexInputFits]->pixels[i].spectro);
							fitsImages[indexInputFits]->pixels[i].spectro = NULL;
							//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
						}							
						free(fitsImages[indexInputFits]->pixels);
						fitsImages[indexInputFits]->pixels = NULL;
					}

					free(resultsInitModelTotal_L[indexInputFits]);		
					free(chisqrfTotal_L[indexInputFits]);
					free(vNumIterTotal_L[indexInputFits]);
					if(configCrontrolFile.SaveSynthesisAdjusted){
						free(vSpectraAjustedTotal_L[indexInputFits]);
					}
					
					printf("\n************************************************************************************************************************\n");
					printf("\n INVERSION OF IMAGE %s ¡¡¡DONE!!!. TIME MAX EXECUTION: %f *********************\n", vInputFileSpectraParalell[indexInputFits].name,vElapsed_execution[indexInputFits]);
					//printf("\n MAX EXECUTION time = %lf seconds\n", vElapsed_execution[indexInputFits]);
					printf("\n TIME TO WRITE FITS IMAGE:  %f seconds to execute \n", timeWriteImage);
					printf("\n************************************************************************************************************************\n\n");
					freeFitsImage(fitsImages[indexInputFits]);
				/*}
				else{
					if(configCrontrolFile.SaveSynthesisAdjusted){
						free(vSpectraAdjustedSplit_L[indexInputFits]);
					}
					free(vSpectraSplit_L[indexInputFits]);
					free(resultsInitModel_L[indexInputFits]);				
					free(vChisqrf_L[indexInputFits]);
					free(vNumIter_L[indexInputFits]);
				}*/
				indexInputFits++;			
			}while(indexInputFits<numFilesPerProcessParallel);
		}

		free(fitsImages);
		free(vMpiRequestScatter);
		free(vMpiRequestInitModel);
		free(vMpiRequestChisqr);
		free(vMpiRequestIter);
		free(vMpiRequestSpectraAd);
		free(vNumPixelsImage);
		free(vMpiRequestReduceExecution);
		free(resultsInitModel_L);
		free(resultsInitModelTotal_L);
		free(chisqrfTotal_L);
		free(vChisqrf_L);
		free(vNumIter_L);
		free(vNumIterTotal_L);
		free(vSpectraSplit_L);
		free(vSpectraAdjustedSplit_L);
		free(vSpectraAjustedTotal_L);
	}


	if(numFilesPerProcess>=1){ // ONE IMAGE PER PROCESSOR
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
			PRECISION timeReadImage	;
			clock_t t, timeTotal;
			t = clock();
			timeTotal = clock();
			fitsImage = readFitsSpectroImage(vInputFileSpectraLocal[indexInputFits].name,0);
			t = clock() - t;
			timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			
			printf("\n************************************************************************************************************************");
			printf("\n\n IDPROC: %d -->  TIME TO READ FITS IMAGE (%s):  %f seconds to execute. NUMBERS OF PIXELS READ %d \n",idProc, vInputFileSpectraLocal[indexInputFits].name, timeReadImage,fitsImage->numPixels); 
			printf("\n************************************************************************************************************************");

			if(fitsImage!=NULL){

				// FITS IMAGE TO STORE OUTPUT PROFILES ADJUSTED
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
						imageStokesAdjust->pixels[i].spectro = calloc (nlambda*NPARMS,sizeof(float));
					}
				}

				//***************************************** INIT MEMORY WITH SIZE OF LAMBDA ****************************************************//
				
				int indexPixel = 0;

				// ALLOCATE MEMORY FOR STORE THE RESULTS 

				vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
				vChisqrf = calloc (fitsImage->numPixels , sizeof(float));
				vNumIter = calloc (fitsImage->numPixels, sizeof(int));
				//t = clock();
				
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
					
					// CLASSICAL ESTIMATES TO GET B, GAMMA, vlos, azimuth
					estimacionesClasicas(wlines[1], vGlobalLambda, nlambda, fitsImage->pixels[indexPixel].spectro, &initModel,1);
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


				clock_t t_write = clock();
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
				timeTotal = clock() - timeTotal;
				t_write = clock() - t_write;
				PRECISION timeTotalExecution = ((PRECISION)timeTotal)/CLOCKS_PER_SEC; // in seconds 
				
				printf("\n************************************************************************************************************************\n");
				printf(" \n IDPROC: %d --> IMAGE INVERSION FOR IMAGE %s ¡¡¡DONE!!!. TIME: %f *********************\n", idProc, vInputFileSpectraLocal[indexInputFits].name,timeTotalExecution);
				printf("\n************************************************************************************************************************\n");				
				
				PRECISION timeToWriteImage = ((PRECISION)t_write)/CLOCKS_PER_SEC; // in seconds 
				
				printf("\n************************************************************************************************************************\n");
				printf(" \n IDPROC: %d --> TIME TO WRITE OUTPUT FILES %s . TIME: %f *********************\n", idProc, vInputFileSpectraLocal[indexInputFits].name,timeToWriteImage);
				printf("\n************************************************************************************************************************\n");				
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


			}
			else{
				printf("\n\n IDPROC: %d --> FITS FILE: %s WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n",idProc, vInputFileSpectraLocal[indexInputFits].name);
				//slog_error(0,"\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
			}

			freeFitsImage(fitsImage);
			//FreeMemoryDerivedSynthesis();
		}
		
		free(vInputFileSpectraLocal);
		free(vOutputNameModelsLocal);
		free(vOutputNameSynthesisAdjustedLocal);
		vInputFileSpectraLocal = NULL;
		vOutputNameModelsLocal = NULL;
		vOutputNameSynthesisAdjustedLocal = NULL;
	}


	if(numFilesPer2ProcessParallel>0){ // CASE DIVIDE EACH IMAGE BETWEEN TWO PROCESS

		// Get the group or processes of the default communicator

		// EACH GROUP PROCESS ONE IMAGE
		MPI_Barrier(vCommunicators[myGroup]);
		numPixels=0;
		clock_t timeTotal;

		if(myGroupRank==groupRoot && access(vInputFileSpectraDiv2Parallel[myGroup].name,F_OK)!=-1){
			
			clock_t t = clock();
			timeTotal = clock();
			fitsImage = readFitsSpectroImage(vInputFileSpectraDiv2Parallel[myGroup].name,1);
			t = clock() - t;
			PRECISION timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			printf("\n TIME TO READ FITS IMAGE %s:  %f seconds to execute . NUMBER OF PIXELS READ: %d \n",vInputFileSpectraDiv2Parallel[myGroup].name, timeReadImage,fitsImage->numPixels); 
			numPixels = fitsImage->numPixels;
		}

		MPI_Barrier(vCommunicators[myGroup]); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY by my group 
		MPI_Bcast(&numPixels, 1, MPI_INT, groupRoot , vCommunicators[myGroup]);
		MPI_Barrier(vCommunicators[myGroup]);



		if(numPixels > 0){
			if(myGroupRank==groupRoot){
				printf("\n***********************  DOING INVERSION: %s *******************************\n\n",vInputFileSpectraDiv2Parallel[myGroup].name );
				resultsInitModelTotal = calloc (numPixels , sizeof(Init_Model));
				chisqrfTotal = calloc (numPixels , sizeof(float));
				vNumIterTotal = calloc (numPixels, sizeof(int));
				if(configCrontrolFile.SaveSynthesisAdjusted)
					vSpectraAjustedTotal = calloc (numPixels*nlambda*NPARMS,sizeof(float));
			}				
			//AllocateMemoryDerivedSynthesis(nlambda);
			int numPixelsProceso = numPixels/myGroupSize;
			int resto = numPixels % myGroupSize;
			int sum = 0;                // Sum of counts. Used to calculate displacements
			int sumSpectro = 0;
			int sumLambda = 0;
			int sendcountsDiv2Pixels [myGroupSize] ; // array describing how many elements to send to each process
			int sendcountsDiv2Spectro [myGroupSize];
			int sendcountsDiv2Lambda [myGroupSize];
			int displsDiv2Pixels [myGroupSize];
			int displsDiv2Spectro [myGroupSize];
			for ( i = 0; i < myGroupSize; i++) {
				sendcountsDiv2Pixels[i] = numPixelsProceso;
				if (resto > 0) {
						sendcountsDiv2Pixels[i]++;
						resto--;
				}
				sendcountsDiv2Spectro[i] = sendcountsDiv2Pixels[i]*nlambda*NPARMS;
				sendcountsDiv2Lambda[i] = sendcountsDiv2Pixels[i]*nlambda;
				displsDiv2Pixels[i] = sum;
				displsDiv2Spectro[i] = sumSpectro;
				//displsLambda[i] = sumLambda;
				sum += sendcountsDiv2Pixels[i];
				sumSpectro += sendcountsDiv2Spectro[i];
				sumLambda += sendcountsDiv2Lambda[i];
			}
			MPI_Barrier(vCommunicators[myGroup]);	
			// SCATTER VPIXELS 
			local_start = MPI_Wtime();
			local_start_scatter = MPI_Wtime();
			vSpectraSplit = calloc(sendcountsDiv2Spectro[myGroupRank],sizeof(float));
			
			if(configCrontrolFile.SaveSynthesisAdjusted)
				vSpectraAdjustedSplit = calloc(sendcountsDiv2Spectro[myGroupRank],sizeof(float));

			if(myGroupRank==groupRoot){
				MPI_Scatterv(fitsImage->spectroImagen, sendcountsDiv2Spectro, displsDiv2Spectro, MPI_FLOAT, vSpectraSplit, sendcountsDiv2Spectro[myGroupRank], MPI_FLOAT, groupRoot, vCommunicators[myGroup]);
			}
			else{
				MPI_Scatterv(NULL, NULL,NULL, MPI_FLOAT, vSpectraSplit, sendcountsDiv2Spectro[myGroupRank], MPI_FLOAT, groupRoot, vCommunicators[myGroup]);
			}	
			local_finish_scatter = MPI_Wtime();

			resultsInitModel = calloc(sendcountsDiv2Pixels[myGroupRank], sizeof(Init_Model));
			vChisqrf = calloc(sendcountsDiv2Pixels[myGroupRank], sizeof(float));
			vNumIter = calloc(sendcountsDiv2Pixels[myGroupRank], sizeof(int));

			local_start_execution = MPI_Wtime();
			for(indexPixel = 0; indexPixel < sendcountsDiv2Pixels[myGroupRank]; indexPixel++){
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

				// CLASSICAL ESTIMATES TO GET B, GAMMA, vlos, azimuth
				estimacionesClasicas(wlines[1], vGlobalLambda, nlambda, vSpectraSplit+(indexPixel*(nlambda*NPARMS)), &initModel,1);

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
			MPI_Gatherv(resultsInitModel, sendcountsDiv2Pixels[myGroupRank], mpiInitModel, resultsInitModelTotal, sendcountsDiv2Pixels, displsDiv2Pixels, mpiInitModel, groupRoot, vCommunicators[myGroup]);
			MPI_Gatherv(vChisqrf, sendcountsDiv2Pixels[myGroupRank], MPI_FLOAT, chisqrfTotal, sendcountsDiv2Pixels, displsDiv2Pixels, MPI_FLOAT, groupRoot, vCommunicators[myGroup]);		
			MPI_Gatherv(vNumIter, sendcountsDiv2Pixels[myGroupRank], MPI_INT, vNumIterTotal, sendcountsDiv2Pixels, displsDiv2Pixels, MPI_INT, groupRoot, vCommunicators[myGroup]);		

			if(configCrontrolFile.SaveSynthesisAdjusted)
				MPI_Gatherv(vSpectraAdjustedSplit, sendcountsDiv2Spectro[myGroupRank], MPI_FLOAT, vSpectraAjustedTotal, sendcountsDiv2Spectro, displsDiv2Spectro, MPI_FLOAT, groupRoot, vCommunicators[myGroup]);

			local_finish_gather = MPI_Wtime();
			local_finish = MPI_Wtime();
			local_elapsed = local_finish - local_start;
			local_elapsed_execution = local_finish_execution - local_start_execution;
			local_elapsed_scatter = local_finish_scatter - local_start_scatter;
			local_elapsed_gather = local_finish_gather - local_start_gather;
			MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, groupRoot, vCommunicators[myGroup]);
			MPI_Reduce(&local_elapsed_execution, &elapsed_execution, 1, MPI_DOUBLE, MPI_MAX, groupRoot, vCommunicators[myGroup]);
			MPI_Reduce(&local_elapsed_scatter, &elapsed_scatter, 1, MPI_DOUBLE, MPI_MAX, groupRoot, vCommunicators[myGroup]);
			MPI_Reduce(&local_elapsed_gather, &elapsed_gather, 1, MPI_DOUBLE, MPI_MAX, groupRoot, vCommunicators[myGroup]);

			if(myGroupRank==groupRoot){
				/*printf("\n Elapsed SCATTER time = %lf seconds\n", elapsed_scatter);
				printf("\n***********************************\n");
				printf("\n Elapsed GATHER time = %lf seconds\n", elapsed_gather);
				printf("\n***********************************\n");			*/
				printf("\n MAX EXECUTION time = %lf seconds\n", elapsed_execution);
				printf("\n***********************************\n");
				/*printf("\n Elapsed TOTAL time = %lf seconds\n", elapsed);
				printf("\n***********************************\n");*/
				double timeWriteImage;
				clock_t t;
				t = clock();	
				if(!writeFitsImageModels(vOutputNameModelsDiv2Parallel[myGroup].name,fitsImage->rows,fitsImage->cols,resultsInitModelTotal,chisqrfTotal,vNumIterTotal,configCrontrolFile.saveChisqr)){
						printf("\n ERROR WRITING FILE OF MODELS: %s",vOutputNameModelsDiv2Parallel[myGroup].name);
				}
				t = clock() - t;
				timeWriteImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n TIME TO WRITE FITS IMAGE:  %f seconds to execute \n", timeWriteImage);
				if(configCrontrolFile.SaveSynthesisAdjusted){
					fitsImage->pixels = calloc(fitsImage->numPixels, sizeof(vpixels));
					for( i=0;i<fitsImage->numPixels;i++){
						fitsImage->pixels[i].spectro = calloc ((fitsImage->numStokes*fitsImage->nLambdas),sizeof(float));
						//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
					}		
					for(indexPixel=0;indexPixel<numPixels;indexPixel++)
					{	
						int kk;
						for (kk = 0; kk < (nlambda * NPARMS); kk++)
						{
							fitsImage->pixels[indexPixel].spectro[kk] = vSpectraAjustedTotal[kk+(indexPixel*(nlambda * NPARMS))] ;
						}
					}					
					// WRITE SINTHETIC PROFILES TO FITS FILE
					if(!writeFitsImageProfiles(vOutputNameSynthesisAdjustedDiv2Parallel[myGroup].name,vInputFileSpectraDiv2Parallel[myGroup].name,fitsImage)){
						printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",vOutputNameSynthesisAdjustedDiv2Parallel[myGroup].name);
					}
					
					for( i=0;i<fitsImage->numPixels;i++){
						free(fitsImage->pixels[i].spectro);
						fitsImage->pixels[i].spectro = NULL;
						//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
					}							
					free(fitsImage->pixels);
					fitsImage->pixels = NULL;
				}
				timeTotal = clock() - timeTotal;
				PRECISION timeTotalExecution = ((PRECISION)timeTotal)/CLOCKS_PER_SEC; // in seconds 
				free(resultsInitModelTotal);		
				free(chisqrfTotal);
				free(vNumIterTotal);
				if(configCrontrolFile.SaveSynthesisAdjusted){
					free(vSpectraAjustedTotal);
				}
				printf("\n************************************************************************************************************************\n");
				printf("\n INVERSION OF IMAGE %s ¡¡¡DONE!!!. TIME: %f *********************\n", vInputFileSpectraDiv2Parallel[myGroup].name,timeTotalExecution);
				printf("\n************************************************************************************************************************\n");					
			}
			else{
				if(configCrontrolFile.SaveSynthesisAdjusted)
					free(vSpectraAdjustedSplit);
				free(vSpectraSplit);
				free(resultsInitModel);				
				free(vChisqrf);
				free(vNumIter);
			}
		}
		else
		{
			if(myGroupRank==groupRoot){
				printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT %s ******************************",vInputFileSpectraDiv2Parallel[myGroup].name);
			}
		}
		if(myGroupRank==groupRoot){
			freeFitsImage(fitsImage);
		}
	}



	FreeMemoryDerivedSynthesis();

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

	if(vInputFileSpectra != NULL) free(vInputFileSpectra);
	if(vInputFileSpectraParalell != NULL) free(vInputFileSpectraParalell);
	if(vInputFileSpectraDiv2Parallel != NULL) free(vInputFileSpectraDiv2Parallel);
	if(vOutputNameModels != NULL) free(vOutputNameModels);
	if(vOutputNameModelsParalell != NULL) free(vOutputNameModelsParalell);
	if(vOutputNameModelsDiv2Parallel != NULL) free(vOutputNameModelsDiv2Parallel);
	if(vOutputNameSynthesisAdjusted != NULL) free(vOutputNameSynthesisAdjusted);
	if(vOutputNameSynthesisAdjustedParallel != NULL) free(vOutputNameSynthesisAdjustedParallel);
	if(vOutputNameSynthesisAdjustedDiv2Parallel != NULL) free(vOutputNameSynthesisAdjustedDiv2Parallel);
	if(vInputFileSpectraLocal != NULL) free(vInputFileSpectraLocal);
	if(vOutputNameModelsLocal != NULL) free(vOutputNameModelsLocal);
	if(vOutputNameSynthesisAdjustedLocal != NULL) free(vOutputNameSynthesisAdjustedLocal);
	
	free(cuantic);
	free(wlines);
	free(vGlobalLambda);
	
	// FREE TYPE OF MPI
	MPI_Type_free(&mpiInitModel);
	MPI_Finalize() ;
	free(G);
	gsl_eigen_symmv_free (workspace);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	return 0;
}

