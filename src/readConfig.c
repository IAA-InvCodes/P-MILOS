#include "readConfig.h"
#include "defines.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <locale.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
//#include <libconfig.h>



/*int readConfigControl(char * configFile, ConfigControl * trolConfig, int printLog){


	config_t cfg;
	config_init(&cfg);
	  
	if(! config_read_file(&cfg, configFile))
	{
		if(printLog) printf("%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	//  number of iterations 
	
	if(config_lookup_int(&cfg, NUMBER_OF_CYCLES, &trolConfig->NumberOfCycles)){
    	if(printLog) printf("NumberOfCycles to apply: %d\n", trolConfig->NumberOfCycles);
	}
  	else
    	if(printLog) printf("No 'NumberOfIterations' setting in configuration file.\n");

	//  file with spectro observed 
	if(config_lookup_string(&cfg, OBSERVED_PROFILES, &trolConfig->ObservedProfiles)){
    	if(printLog) printf("ObservedProfiles to apply: %s\n", trolConfig->ObservedProfiles);
	}
  	else
    	if(printLog) printf("No 'ObservedProfiles' setting in configuration file.\n");

	//  file with stray light 
	if(config_lookup_string(&cfg, STRAY_LIGHT_FILE, &trolConfig->StrayLightFile)){
    	if(printLog) printf("StrayLightFile to apply: %s\n", trolConfig->StrayLightFile);
	}
  	else
    	if(printLog) printf("No 'StrayLightFile' setting in configuration file.\n");

	//  file with psf of instrument 
	if(config_lookup_string(&cfg, PSF_FILE, &trolConfig->PSFFile)){
    	if(printLog) printf("PSFFile to apply: %s\n", trolConfig->PSFFile);
	}
  	else
    	if(printLog) printf("No 'PSFFile' setting in configuration file.\n");		

	//  fits file with wave lenght 
	if(config_lookup_string(&cfg, WAVE_LENGHT_FILE, &trolConfig->WavelengthFile)){
    	if(printLog) printf("WavelengthFile to apply: %s\n", trolConfig->WavelengthFile);
	}
  	else
    	if(printLog) printf("No 'WavelengthFile' setting in configuration file.\n");		

	//  file with atomic lines 
	if(config_lookup_string(&cfg, ATOMIC_PARAMETERS_FILE, &trolConfig->AtomicParametersFile)){
    	if(printLog) printf("AtomicParametersFile to apply: %s\n", trolConfig->AtomicParametersFile);
	}
  	else
    	if(printLog) printf("No 'AtomicParametersFile' setting in configuration file.\n");

	//  READ INITIAL GUESS MODEL 
	if(config_lookup_string(&cfg, INITIAL_GUESS_MODEL, &trolConfig->InitialGuessModel)){
    	if(printLog) printf("InitialGuessModel to apply: %s\n", trolConfig->InitialGuessModel);
	}
  	else
    	if(printLog) printf("No 'InitialGuessModel' setting in configuration file.\n");		

	//  Weight for stokes I  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESI, &trolConfig->WeightForStokes[0])){
    	if(printLog) printf("WeightForStokesI to apply: %f\n", trolConfig->WeightForStokes[0]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesI' setting in configuration file.\n");

	//  Weight for stokes Q  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESQ, &trolConfig->WeightForStokes[1])){
    	if(printLog) printf("WeightForStokesQ to apply: %f\n", trolConfig->WeightForStokes[1]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesQ' setting in configuration file.\n");		

	//  Weight for stokes U  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESU, &trolConfig->WeightForStokes[2])){
    	if(printLog) printf("WeightForStokesU to apply: %f\n", trolConfig->WeightForStokes[2]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesU' setting in configuration file.\n");		

	//  Weight for stokes V 
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESV, &trolConfig->WeightForStokes[3])){
    	if(printLog) printf("WeightForStokesV to apply: %f\n", trolConfig->WeightForStokes[3]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesV' setting in configuration file.\n");

	//  invert macroturbulence
	if(config_lookup_bool(&cfg, INVERT_MACROTURBULENCE, &trolConfig->InvertMacroturbulence)){
    	if(printLog) printf("InvertMacroturbulence to apply: %d\n", trolConfig->InvertMacroturbulence);
	}
  	else
    	if(printLog) printf("No 'InvertMacroturbulence' setting in configuration file.\n");		

	//  invert filling factor 
	if(config_lookup_bool(&cfg, INVERT_FILLING_FACTOR, &trolConfig->InvertFillingFactor)){
    	if(printLog) printf("InvertFillingFactor to apply: %d\n", trolConfig->InvertFillingFactor);
	}
  	else
    	if(printLog) printf("No 'InvertFillingFactor' setting in configuration file.\n");		

	//  invert stray light factor
	if(config_lookup_bool(&cfg, INVERT_STRAY_LIGHT_FACTOR, &trolConfig->InvertStrayLightFactor)){
    	if(printLog) printf("InvertStrayLightFactor to apply: %d\n", trolConfig->InvertStrayLightFactor);
	}
  	else
    	if(printLog) printf("No 'InvertStrayLightFactor' setting in configuration file.\n");

	//  Read mu : Scalar containing the cosine of the heliocentric angle 
	if(config_lookup_float(&cfg, MU, &trolConfig->mu)){
    	if(printLog) printf("mu to apply: %f\n", trolConfig->mu);
	}
  	else
    	if(printLog) printf("No 'mu' setting in configuration file.\n");		

	//  Estimated S/N FOR I 
	if(config_lookup_int(&cfg, ESTIMATEDSNFORI, &trolConfig->EstimatedSNForI)){
    	if(printLog) printf("EstimatedSNForI to apply: %d\n", trolConfig->EstimatedSNForI);
	}
  	else
    	if(printLog) printf("No 'EstimatedSNForI' setting in configuration file.\n");		

	//  Continuum Contrast
	if(config_lookup_int(&cfg, CONTINUUM_CONTRAST, &trolConfig->ContinuumContrast)){
    	if(printLog) printf("ContinuumContrast to apply: %d\n", trolConfig->ContinuumContrast);
	}
  	else
    	if(printLog) printf("No 'ContinuumContrast' setting in configuration file.\n");			

	//  initial diagonal element
	if(config_lookup_float(&cfg, INITIAL_DIAGONAL_ELEMENT, &trolConfig->InitialDiagonalElement)){
    	if(printLog) printf("InitialDiagonalElement to apply: %f\n", trolConfig->InitialDiagonalElement);
	}
  	else
    	if(printLog) printf("No 'InitialDiagonalElement' setting in configuration file.\n");	
			

	//  Use interpolar with splines or linear
	if(config_lookup_int(&cfg, USE_INTERPOLAR_SPLINES_OR_LINEAR, &trolConfig->useInterpolarSplinesOrLinear)){
    	if(printLog) printf("useInterpolarSplinesOrLinear to apply: %d\n", trolConfig->useInterpolarSplinesOrLinear);
	}
  	else
    	if(printLog) printf("No 'useInterpolarSplinesOrLinear' setting in configuration file.\n");					

	//  Indicate if convolve with PSF
	if(config_lookup_bool(&cfg, CONVOLVE_WITH_PSF, &trolConfig->ConvolveWithPSF)){
    	if(printLog) printf("ConvolveWithPSF to apply: %d\n", trolConfig->ConvolveWithPSF);
	}
  	else
    	if(printLog) printf("No 'ConvolveWithPSF' setting in configuration file.\n");							

	//  FWHM
	if(config_lookup_float(&cfg, FWHM_FILE, &trolConfig->FWHM)){
    	if(printLog) printf("FWHM to apply: %f\n", trolConfig->FWHM);
	}
  	else
    	if(printLog) printf("No 'FWHM' setting in configuration file.\n");			

	//  TYPE CONVOLUTION TO USE 
	if(config_lookup_string(&cfg, TYPE_CONVOLUTION, &trolConfig->TypeConvolution)){
    	if(printLog) printf("TypeConvolution to apply: %s\n", trolConfig->TypeConvolution);
	}
  	else
    	if(printLog) printf("No 'TypeConvolution' setting in configuration file.\n");			

	//  Gas pressure at surface 1
	if(config_lookup_float(&cfg, GAS_PRESSURE_AT_SURFACE_1, &trolConfig->GasPressureAtSurface1)){
    	if(printLog) printf("GasPressureAtSurface1to apply: %f\n", trolConfig->GasPressureAtSurface1);
	}
  	else
    	if(printLog) printf("No 'GasPressureAtSurface1' setting in configuration file.\n");					

	//  Gas pressure at surface 2
	if(config_lookup_float(&cfg, GAS_PRESSURE_AT_SURFACE_2, &trolConfig->GasPressureAtSurface2)){
    	if(printLog) printf("GasPressureAtSurface2 to apply: %f\n", trolConfig->GasPressureAtSurface2);
	}
  	else
    	if(printLog) printf("No 'GasPressureAtSurface2' setting in configuration file.\n");			

	//  Magnetic pressure term
	if(config_lookup_float(&cfg, MAGNETIC_PRESSURE_TERM, &trolConfig->MagneticPressureTerm)){
    	if(printLog) printf("MagneticPressureTerm to apply: %f\n", trolConfig->MagneticPressureTerm);
	}
  	else
    	if(printLog) printf("No 'MagneticPressureTerm' setting in configuration file.\n");		

	//  Numbe of lines tu use
	if(config_lookup_int(&cfg, NTL, &trolConfig->ntl)){
    	if(printLog) printf("ntl to apply: %d\n", trolConfig->ntl);
	}
  	else
    	if(printLog) printf("No 'ntl' setting in configuration file.\n");		

	//  Number of wavelenght observed
	if(config_lookup_int(&cfg, NLIOBS, &trolConfig->nliobs)){
    	if(printLog) printf("nliobs to apply: %d\n", trolConfig->nliobs);
	}
  	else
    	if(printLog) printf("No 'nliobs' setting in configuration file.\n");	

	//  Central wavelenght 
	if(config_lookup_float(&cfg, CENTRAL_WAVE_LENGHT, &trolConfig->CentralWaveLenght)){
    	if(printLog) printf("CentralWaveLenght to apply: %f\n", trolConfig->CentralWaveLenght);
	}
  	else
    	if(printLog) printf("No 'CentralWaveLenght' setting in configuration file.\n");	

	//  indicate if calculate eta0 for compute output model
	if(config_lookup_int(&cfg, ETA0_LINE_TO_CONTINUUM_ABSORPTION, &trolConfig->fix[0])){
    	if(printLog) printf("ETA0 to apply: %d\n", trolConfig->fix[0]);
	}
  	else
    	if(printLog) printf("No 'ETA0' setting in configuration file.\n");	

	//  indicate if calculate B for compute output model
	if(config_lookup_int(&cfg, B_MAGNETIC_FIELD_STRENGTH, &trolConfig->fix[1])){
    	if(printLog) printf("B to apply: %d\n", trolConfig->fix[1]);
	}
  	else
    	if(printLog) printf("No 'B' setting in configuration file.\n");	

	//  Indicate if calculate VLOS for output model
	if(config_lookup_int(&cfg, VLOS_LINE_OF_SIGHT_VELOCITY, &trolConfig->fix[2])){
    	if(printLog) printf("VLOS to apply: %d\n", trolConfig->fix[2]);
	}
  	else
    	if(printLog) printf("No 'VLOS' setting in configuration file.\n");			

	//  Indicate if calculate DOPP 
	if(config_lookup_int(&cfg, DOPP_DOOPLER_WIDTH, &trolConfig->fix[3])){
    	if(printLog) printf("DOPP to apply: %d\n", trolConfig->fix[3]);
	}
  	else
    	if(printLog) printf("No 'DOPP' setting in configuration file.\n");		

	//  indicate if calculate aa  
	if(config_lookup_int(&cfg, AA_DAMPING_PARAMETER, &trolConfig->fix[4])){
    	if(printLog) printf("AA to apply: %d\n", trolConfig->fix[4]);
	}
  	else
    	if(printLog) printf("No 'AA' setting in configuration file.\n");		

	//  indicate if calculate gm
	if(config_lookup_int(&cfg, GM_MAGNETIC_FIELD_INCLINATION, &trolConfig->fix[5])){
    	if(printLog) printf("GM to apply: %d\n", trolConfig->fix[5]);
	}
  	else
    	if(printLog) printf("No 'GM' setting in configuration file.\n");		

	//  indicate if calculate azimuth 
	if(config_lookup_int(&cfg, AZ_MAGNETIC_FIELD_AZIMUTH, &trolConfig->fix[6])){
    	if(printLog) printf("AZ to apply: %d\n", trolConfig->fix[6]);
	}
  	else
    	if(printLog) printf("No 'AZ' setting in configuration file.\n");		

	//  indicate if calculate s0
	if(config_lookup_int(&cfg, S0_SOURCE_FUNCTION_CONSTANT, &trolConfig->fix[7])){
    	if(printLog) printf("S0 to apply: %d\n", trolConfig->fix[7]);
	}
  	else
    	if(printLog) printf("No 'S0' setting in configuration file.\n");	

	//  indicate if calculate s1 
	if(config_lookup_int(&cfg, S1_SOURCE_FUNCTION_GRADIENT, &trolConfig->fix[8])){
    	if(printLog) printf("S1 to apply: %d\n", trolConfig->fix[8]);
	}
  	else
    	if(printLog) printf("No 'S1' setting in configuration file.\n");						

	//  indicate if calculate macroturbulence
	if(config_lookup_int(&cfg, MAC_MACROTURBULENT_VELOCITY, &trolConfig->fix[9])){
    	if(printLog) printf("MAC to apply: %d\n", trolConfig->fix[9]);
	}
  	else
    	if(printLog) printf("No 'MAC' setting in configuration file.\n");			

	//  indicate if calculate filling factor 
	if(config_lookup_int(&cfg, ALPHA_FILLING_FACTOR, &trolConfig->fix[10])){
    	if(printLog) printf("ALPHA to apply: %d\n", trolConfig->fix[10]);
	}
  	else
    	if(printLog) printf("No 'ALPHA' setting in configuration file.\n");	

	// if save chisqr 
	if(config_lookup_bool(&cfg, SAVE_CHISQR, &trolConfig->saveChisqr)){
    	if(printLog) printf("saveChisqr to apply: %d\n", trolConfig->saveChisqr);
	}
  	else{
    	if(printLog) printf("No 'saveChisqr' setting in configuration file. Default = FALSE\n");	
		trolConfig->saveChisqr=0;
	}

	//  if use classical estimates
	if(config_lookup_bool(&cfg, USE_CLASSICAL_ESTIMATES, &trolConfig->UseClassicalEstimates)){
    	if(printLog) printf("UseClassicalEstimates to apply: %d\n", trolConfig->UseClassicalEstimates);
	}
  	else
    	if(printLog) printf("No 'UseClassicalEstimates' setting in configuration file.\n");	

	//  if use rte inversion 
	if(config_lookup_bool(&cfg, USE_RTE_INVERSION, &trolConfig->UseRTEInversion)){
    	if(printLog) printf("UseRTEInversion to apply: %d\n", trolConfig->UseRTEInversion);
	}
  	else
    	if(printLog) printf("No 'UseRTEInversion' setting in configuration file.\n");	

	//  if save synthesis profile
	if(config_lookup_bool(&cfg, SAVE_SYNTHESIS_PROFILE, &trolConfig->SaveSynthesisAdjusted)){
    	if(printLog) printf("SaveSynthesisAdjusted to apply: %d\n", trolConfig->SaveSynthesisAdjusted);
	}
  	else
    	if(printLog) printf("No 'SaveSynthesisAdjusted' setting in configuration file.\n");	

	//  Output model file 
	if(config_lookup_string(&cfg, OUTPUT_MODEL_FILE, &trolConfig->OutputModelFile)){
    	if(printLog) printf("OutputModelFile to apply: %s\n", trolConfig->OutputModelFile);
	}
  	else
    	if(printLog) printf("No 'OutputModelFile' setting in configuration file.\n");	

	//  Output synthesis file 
	if(config_lookup_string(&cfg, OUTPUT_SYNTHESIS_FILE, &trolConfig->OutputSynthesisFile)){
    	if(printLog) printf("OutputSynthesisFile to apply: %s\n", trolConfig->OutputSynthesisFile);
	}
  	else
    	if(printLog) printf("No 'OutputSynthesisFile' setting in configuration file.\n");	

	//  noise to apply in all stoke parameter, only will be applied if sigma is empty
	if(config_lookup_float(&cfg, NOISE_FILE, &trolConfig->noise)){
    	if(printLog) printf("noise to apply: %f\n", trolConfig->noise);
	}
  	else
    	if(printLog) printf("No 'noise' setting in configuration file.\n");

	//  sigma to apply in each stokes parameter 
	if(config_lookup_float(&cfg, SIGMA_FILE, trolConfig->sigma)){
    	if(printLog) printf("Sigma to apply\n");
	}
  	else{
   	if(printLog) printf("No 'Sigma' setting in configuration file. Assigning noise to each parameter of the vector. \n");			
		trolConfig->sigma[0] = trolConfig->noise;
		trolConfig->sigma[1] = trolConfig->noise;
		trolConfig->sigma[2] = trolConfig->noise;
		trolConfig->sigma[3] = trolConfig->noise;
	}

		

	//  Optional minimum relative difference between two succesive merit-function values
	if(config_lookup_float(&cfg, TOPLIM_FILE, &trolConfig->toplim)){
    	if(printLog) printf("toplim to apply: %e\n", trolConfig->toplim);
	}
  	else
    	if(printLog) printf("No 'toplim' setting in configuration file.\n");					

	return 1;

}*/

/**
 * Read Cuantic data from a file with the Cuantic Lines. 
 * @param inputLineFile
 * @param cuanticDat
 * @param line2Read
 * @return Returns the value of central wavelength if the line is found or 0 in other case. 
 * 
 * */
PRECISION readFileCuanticLines(const char * inputLineFile, PRECISION * cuanticDat, int line2Read, int printLog){
	// try open the file with the 
	FILE * fp;
	char * line = NULL;

   size_t len = 0;
   ssize_t read;
	char atomo [2];
	fp = fopen(inputLineFile, "r");
   if(fp == NULL)	return 0;

	int indexLine, ionicState;
	int found = 0;
	double damping, potentialExcitation, logGf;
	PRECISION lambdaLine;
	int SLOI, SUPI;
	PRECISION LLOI,JLOI,LUPI,JUPI;
	char levelD,levelU;

	while ((read = getline(&line, &len, fp)) != -1  && !found ) {
		sscanf(line,"%i=%s %i %lf %lf %lf -%lf %i%c %lf- %i%c %lf",&indexLine,atomo,&ionicState,&lambdaLine,&damping,&potentialExcitation,&logGf,&SLOI,&levelD,&JLOI,&SUPI,&levelU,&JUPI);
		if(indexLine==line2Read){ // read the rest of the line, else read next line
			switch (levelD)
			{
			case 'S':
				LLOI = 0;
				break;
			case 'P':
				LLOI = 1;
				break;
			case 'D':
				LLOI = 2;
				break;				
			case 'F':
				LLOI = 3;
				break;
			case 'G':
				LLOI = 4;
				break;				
			case 'H':
				LLOI = 5;
				break;
			case 'J':
				LLOI = 6;
				break;
			default:
				break;
			}
			switch (levelU)
			{
			case 'S':
				LUPI = 0;
				break;
			case 'P':
				LUPI = 1;
				break;
			case 'D':
				LUPI = 2;
				break;				
			case 'F':
				LUPI = 3;
				break;
			case 'G':
				LUPI = 4;
				break;				
			case 'H':
				LUPI = 5;
				break;
			case 'J':
				LUPI = 6;
				break;
			default:
				break;
			}
			SLOI= (SLOI-1)/2;
			SUPI= (SUPI-1)/2;

			found = 1; 
		}
   }
	cuanticDat[0] =1 ; // LINE NUMBER 1
	cuanticDat[1] = SLOI;
	cuanticDat[2] = LLOI;
	cuanticDat[3] = JLOI;
	cuanticDat[4] = SUPI;
	cuanticDat[5] = LUPI;
	cuanticDat[6] = JUPI;
	
	if(printLog){
		printf("\n***********************************");
		printf("\n\n QUANTUM NUMBERS READ FROM FILE, FOR CENTRAL WAVELENGTH %lf: \n",lambdaLine);
		printf("\n\tSLOI: %fd",cuanticDat[1]);
		printf("\n\tLLOI: %fd",cuanticDat[2]);
		printf("\n\tJLOI: %fd",cuanticDat[3]);
		printf("\n\tSUPI: %fd",cuanticDat[4]);
		printf("\n\tLUPI: %fd",cuanticDat[5]);
		printf("\n\tJUPI: %fd",cuanticDat[6]);
		
	}
	if(!found)
		return 0;
	else
		return lambdaLine;

}


int readInitialModel(Init_Model * INIT_MODEL, char * fileInitModel){
	
	FILE * fReadInitModel;
	char * line = NULL;
	size_t len = 0;
   ssize_t read;
	char comment[200], name[100];
	
	fReadInitModel = fopen(fileInitModel, "r");
	if (fReadInitModel == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileInitModel);
		return 0;
	}
	
	while ((read = getline(&line, &len, fReadInitModel)) != -1) {
		double aux_value;
		sscanf(line,"%99[^:]:%lf%99[^!]!",name, &aux_value,comment);
		if(strstr(name,"INITIAL_MODEL_B")!=NULL){ // B
			INIT_MODEL->B = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_GM")!=NULL){ // GM
			INIT_MODEL->gm = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_AZI")!=NULL){ // AZI
			INIT_MODEL->az = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_ETHA0")!=NULL){ // ETHA0
			INIT_MODEL->eta0 = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_LAMBDADOPP")!=NULL){ // LAMBDADOPP
			INIT_MODEL->dopp = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_AA")!=NULL){ // AA
			INIT_MODEL->aa = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_ALFA")!=NULL){ // ALFA
			INIT_MODEL->alfa = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_MAC")!=NULL){ // MAC
			INIT_MODEL->mac = aux_value;
		}		
		if(strstr(name,"INITIAL_MODEL_VLOS")!=NULL){ // VLOS
			INIT_MODEL->vlos = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_S0")!=NULL){ // S0
			INIT_MODEL->S0 = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_S1")!=NULL){ // S1
			INIT_MODEL->S1 = aux_value;
		}				
	}
	fclose(fReadInitModel);
	return 1;
}

/*int readInitialModel(Init_Model * INIT_MODEL,const char * fileInitModel){
	
	config_t cfg;
  	
	config_init(&cfg);
	if(! config_read_file(&cfg, fileInitModel))
	{
		printf("%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	//  ETA0
	if(!config_lookup_float(&cfg, INITIAL_MODEL_ETHA0, &INIT_MODEL->eta0))
    	printf( "No 'ETA0' setting in configuration file. Used by default \n");

	//  B
	if(!config_lookup_float(&cfg, INITIAL_MODEL_B, &INIT_MODEL->B))
    	printf("No 'B' setting in configuration file. Used by default \n");

	//  vlos
	if(!config_lookup_float(&cfg, INITIAL_MODEL_VLOS, &INIT_MODEL->vlos))
    	printf("No 'VLOS' setting in configuration file. Used by default \n");

	//  LAMBDADOPP
	if(!config_lookup_float(&cfg, INITIAL_MODEL_LAMBDADOPP, &INIT_MODEL->dopp))
    	printf("No 'LAMBDADOPP' setting in configuration file. Used by default \n");

	//  AA
	if(!config_lookup_float(&cfg, INITIAL_MODEL_AA, &INIT_MODEL->aa))
    	printf("No 'AA' setting in configuration file. Used by default \n");		

	//  GM
	if(!config_lookup_float(&cfg, INITIAL_MODEL_GM, &INIT_MODEL->gm))
    	printf("No 'GM' setting in configuration file. Used by default \n");		

	//  AZI
	if(!config_lookup_float(&cfg, INITIAL_MODEL_AZI, &INIT_MODEL->az))
    	printf("No 'AZI' setting in configuration file. Used by default \n");		

	//  S0
	if(!config_lookup_float(&cfg, INITIAL_MODEL_S0, &INIT_MODEL->S0))
    	printf("No 'S0' setting in configuration file. Used by default \n");				

	//  S1
	if(!config_lookup_float(&cfg, INITIAL_MODEL_S1, &INIT_MODEL->S1))
    	printf("No 'S1' setting in configuration file. Used by default \n");				

	//  MAC
	if(!config_lookup_float(&cfg, INITIAL_MODEL_MAC, &INIT_MODEL->mac))
    	printf("No 'MAC' setting in configuration file. Used by default \n");				

	//  ALFA
	if(!config_lookup_float(&cfg, INITIAL_MODEL_ALFA, &INIT_MODEL->alfa))
    	printf("No 'ALFA' setting in configuration file. Used by default \n");						

	config_destroy(&cfg);

	return EXIT_SUCCESS;
}
*/


/**
 * Read Cuantic data from a file with the Cuantic Lines. 
 * @param fileMallaGrid file with malla drid 
 * @param initialLambda variable to store the initial lambda 
 * @param step variable to store each step
 * @param finalLambda variable to store finallambda
 * @param printlog decide if print log or not
 * @return Returns number of the index line to read from cuantic number lines, return 0 if line is not found
 * 
 * 
 * */
int readMallaGrid(const char * fileMallaGrid, PRECISION * initialLambda, PRECISION * step, PRECISION * finalLambda, int printLog){
	// try open the file with the 
	FILE * fp;
	char * line = NULL;

   size_t len = 0;
   ssize_t read;
	
	fp = fopen(fileMallaGrid, "r");
   if(fp == NULL)	return 0;

	int indexLine;
	int found = 0, dataRead = 0;;
	
	
	char name[100];

	while ((read = getline(&line, &len, fp)) != -1 && !dataRead){
		//printf("\n linea leida %s",line);
		if(found){ //1                       :        -624.37,        21.53,     1765.46
			//sscanf(line,"%i,%i%99[^:]:%lf,%lf,%lf",&indexLine,initialLambda,step,finalLambda);
			sscanf(line,"%i%99[^:]:%lf%lf%lf",&indexLine,name,initialLambda,step,finalLambda);
			dataRead = 1;
		}
		else{
			if(strncmp(line,"------",6)==0){
				found = 1;
			}
		}
	}

	if(printLog){	
		printf("\n\n MALLA GRID --- indexline %d , initial lambda: %lf, step: %lf, finallambda: %lf\n",indexLine,*initialLambda,*step,*finalLambda);
	}

	if(dataRead)
		return indexLine;
	else
		return 0;
	

}


/**
 * 
 * */
int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, const char * nameInputPSF, PRECISION centralWaveLenght){

	// first of all read the number of lines to give size for arrays deltaLambda and PSF
	FILE *fp;

	// alloc memory 

	char * line = NULL;
	size_t len = 0;
   ssize_t read;
	fp=fopen(nameInputPSF,"r");
	if(fp==NULL)
	{
	printf("File \"%s\" does not exist!!!\n",nameInputPSF);
			return 0;
	}	
	int index =0;
	while ((read = getline(&line, &len, fp)) != -1) {
		double delta, psf;
		sscanf(line,"%le  %le", &delta, &psf);
		//deltaLambda[index] = (delta/1000)+centralWaveLenght;
		deltaLambda[index] = delta;
		PSF[index] = psf;
		index++;
	}

	fclose(fp);
	return 1;
}


void loadInitialValues(ConfigControl * configControlFile){

	// array of weight 
	configControlFile->WeightForStokes[0]=1.;
	configControlFile->WeightForStokes[1]=1.;
	configControlFile->WeightForStokes[2]=1.;
	configControlFile->WeightForStokes[3]=1.;

	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	int i;
	for(i=0;i<11;i++){
		configControlFile->fix[i]= 0;
		configControlFile->fix2[i]= 0;
	}
	/*configControlFile->fix[10] = 0;
	configControlFile->fix2[10]= 0;*/

	configControlFile->noise = NOISE_SIGMA;
	configControlFile->sigma[0] = NOISE_SIGMA;
	configControlFile->sigma[1] = NOISE_SIGMA;
	configControlFile->sigma[2] = NOISE_SIGMA;
	configControlFile->sigma[3] = NOISE_SIGMA;

	
	configControlFile->InitialDiagonalElement = ILAMBDA;
	configControlFile->toplim = TOPLIM;
	configControlFile->mu = AH;
	configControlFile->saveChisqr = 1;
	configControlFile->SaveSynthesisAdjusted=1;

	configControlFile->subx1 = 0;
	configControlFile->subx2 = 0;
	configControlFile->suby1 = 0;
	configControlFile->suby2 = 0;

	configControlFile->useFFT = 0; // by default direct convolution

}

int readParametersFileInput(char * fileParameters,  ConfigControl * trolConfig, int printLog){

	// try open the file with the 
	FILE * fReadParameters;
	float aux;
	char LINE [4096], * returnLine;
	char comment[200], name[100];
	fReadParameters = fopen(fileParameters, "r");
	if (fReadParameters == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;
	}
	int rfscanf; 
	
	/***************************  NUMBER OF CYCLES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->NumberOfCycles,comment);
	
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NumberOfCycles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	/*if( trolConfig->NumberOfCycles < 0){
		printf("milos: Error in NumberOfCycles parameter. review it. Not accepted: %d\n", trolConfig->NumberOfCycles);
		return 0;
	}*/
	if(printLog) printf("NumberOfCycles to apply: %i\n", trolConfig->NumberOfCycles);

	/***************************  OBSERVED PROFILES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->ObservedProfiles,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Observed Profiles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Observed profiles to apply: %s\n", trolConfig->ObservedProfiles);

	/***************************  STRAY LIGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->StrayLightFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Stray light file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Stray light file to apply: %s\n", trolConfig->StrayLightFile);


	/***************************  PSF FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->PSFFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param PSF file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("PSF file to apply: %s\n", trolConfig->PSFFile);

	/*************************** WAVELENGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->WavelengthFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param wavelength file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("wavelength file to apply: %s\n", trolConfig->WavelengthFile);

	/*************************** WAVELENGHT  GRID FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->MallaGrid,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param wavelength GRID file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("wavelength file grid to apply: %s\n", trolConfig->MallaGrid);


	/*************************** ATOMIC PARAMETER  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->AtomicParametersFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Atomic parameters file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Atomic parameters file to apply: %s\n", trolConfig->AtomicParametersFile);

	/*************************** INITIAL GUESS MODEL   FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial guess model 1 to apply: %s\n", trolConfig->InitialGuessModel);

	/*************************** INITIAL GUESS MODEL  2  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel_2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial guess model 2 to apply: %s\n", trolConfig->InitialGuessModel_2);

	/*************************** WEIGHT FOT STOKES I ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	trolConfig->WeightForStokes[0] = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes I to apply: %lf\n", trolConfig->WeightForStokes[0]);

	/*************************** WEIGHT FOT STOKES Q ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	trolConfig->WeightForStokes[1] = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes Q. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes Q to apply: %lf\n", trolConfig->WeightForStokes[1]);

	/*************************** WEIGHT FOT STOKES U ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	trolConfig->WeightForStokes[2] = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes U. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes U to apply: %lf\n", trolConfig->WeightForStokes[2]);

	/*************************** WEIGHT FOT STOKES V ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	trolConfig->WeightForStokes[3] = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes V. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes V to apply: %lf\n", trolConfig->WeightForStokes[3]);


	/*************************** MU ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->mu,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param mu=cos (theta). Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("mu=cos (theta) to apply: %f\n", trolConfig->mu);

	/*************************** EstimatedSNForI ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->EstimatedSNForI,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Estimated S/N for I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog){
		printf("Estimated S/N for I  to apply: %i\n", trolConfig->EstimatedSNForI);
		printf("Estimated noise for I  to apply: %lf\n", 1.0/trolConfig->EstimatedSNForI);
	} 
	trolConfig->noise = 1.0/trolConfig->EstimatedSNForI;
	trolConfig->sigma[0] = trolConfig->noise*trolConfig->noise;
	trolConfig->sigma[1] = trolConfig->sigma[0];
	trolConfig->sigma[2] = trolConfig->sigma[0];
	trolConfig->sigma[3] = trolConfig->sigma[0];
	// PUT VALUES IN ARRAY OF SIGMA 


 	/*************************** CONTINIUM CONTRAST ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ContinuumContrast,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Continuum contrast. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Continuum contrast to apply: %i\n", trolConfig->ContinuumContrast);

	/*************************** INITIAL_DIAGONAL_ELEMENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux ,comment);
	trolConfig->InitialDiagonalElement = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial diagonal element. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial diagonal element  to apply: %le\n", trolConfig->InitialDiagonalElement);

	/*************************** USE_INTERPOLAR_SPLINES_OR_LINEAR ********************************************/
	
	/*returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->useInterpolarSplinesOrLinear,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Splines/Linear Interpolation. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Splines/Linear Interpolation to apply: %d\n", trolConfig->useInterpolarSplinesOrLinear);*/

	/*************************** USE PSF FILTER GAUSSIAN ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ConvolveWithPSF,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use Gaussian PSF . Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use Gaussian PSF  to apply: %d\n", trolConfig->ConvolveWithPSF);

	/*************************** FWHM ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->FWHM,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param FWHM. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("FWHM to apply: %f\n", trolConfig->FWHM);

	/*************************** NTL ********************************************/
	
	/*returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ntl,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NTL. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NTL to apply: %d\n", trolConfig->ntl);*/

	/*************************** NLIOBS ********************************************/
	
	/*returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->nliobs,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NLIOBS. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NLIOBS to apply: %d\n", trolConfig->nliobs);*/


	/*************************** OUTFILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->outputPrefix,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param OUTFILE. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("OUTFILE PREFIX to apply: %s\n", trolConfig->outputPrefix);


	/*************************** CENTRAL_WAVE_LENGHT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->CentralWaveLenght,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Central Wave Lenght. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Central Wave Lenght  to apply: %f\n", trolConfig->CentralWaveLenght);

	/*************************** ETA0_LINE_TO_CONTINUUM_ABSORPTION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[0],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Line To Continiuum Absorption. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Line To Continiuum Absorption to apply: %d\n", trolConfig->fix[0]);

	/*************************** B_MAGNETIC_FIELD_STRENGTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[1],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic Field Strength. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic Field Strength to apply: %d\n", trolConfig->fix[1]);

	/*************************** VLOS_LINE_OF_SIGHT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[2],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Line Of Sight Velocity. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Line Of Sight Velocity to apply: %d\n", trolConfig->fix[2]);

	/*************************** DOPP_DOOPLER_WIDTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[3],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Doopler Width. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Doopler Width to apply: %d\n", trolConfig->fix[3]);

	/*************************** AA_DAMPING_PARAMETER ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[4],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Damping Parameter. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Damping Parameter to apply: %d\n", trolConfig->fix[4]);

	/*************************** GM_MAGNETIC_FIELD_INCLINATION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[5],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic FieldInclination. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic FieldInclination to apply: %d\n", trolConfig->fix[5]);

	/*************************** AZ_MAGNETIC_FIELD_AZIMUTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[6],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic Field Azimuth. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic Field Azimuth to apply: %d\n", trolConfig->fix[6]);

	/*************************** S0_SOURCE_FUNCTION_CONSTANT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[7],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Source Function Constant. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Source Function Constant to apply: %d\n", trolConfig->fix[7]);

	/*************************** S1_SOURCE_FUNCTION_GRADIENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[8],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Source Function Gradient . Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Source Function Gradient  to apply: %d\n", trolConfig->fix[8]);

	/*************************** MAC_MACROTURBULENT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[9],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Macroturbulent Velocity. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Macroturbulent Velocity to apply: %d\n", trolConfig->fix[9]);

	/*************************** ALPHA_FILLING_FACTOR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[10],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Filling Factor. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Filling Factor  to apply: %d\n", trolConfig->fix[10]);

	/*************************** SAVE_CHISQR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->saveChisqr,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Chisqr. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Save Chisqr  to apply: %d\n", trolConfig->saveChisqr);

	/*************************** USE_CLASSICAL_ESTIMATES ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->UseClassicalEstimates,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use Classical Estimates. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use Classical Estimates to apply: %d\n", trolConfig->UseClassicalEstimates);

	/*************************** USE_RTE_INVERSION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->UseRTEInversion,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use RTE Inversion. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use RTE Inversion to apply: %d\n", trolConfig->UseRTEInversion);


	/*************************** SAVE SYNTHESIS PROFILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->SaveSynthesisAdjusted,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Synthesis Profile Adjusted. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Save Synthesis Profile to apply: %d\n", trolConfig->SaveSynthesisAdjusted);

	/*************************** OUTPUT_MODEL_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->OutputModelFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Output Model File. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Output Model File  to apply: %s\n", trolConfig->OutputModelFile);

	/*************************** Type Output File Model  ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->typeFileOutputModel,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param type of output file model. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Output File Model to apply: %d\n", trolConfig->typeFileOutputModel);


	/*************************** OUTPUT_SYNTHESIS_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->OutputSynthesisFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Output Synthesis File. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Output Synthesis File to apply: %s\n", trolConfig->OutputSynthesisFile);

	/*************************** SIGMA_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	trolConfig->noise = aux;
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NOISE. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NOISE to apply: %le\n", trolConfig->noise);


	/*************************** TOPLIM_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->toplim,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param TOPLIM. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("TOPLIM to apply: %le\n", trolConfig->toplim);

	/*************************** NX ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->nx,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NX. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NX to apply: %d\n", trolConfig->nx);

	/*************************** NY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ny,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NY. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NY to apply: %d\n", trolConfig->ny);

	/*************************** subx1 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->subx1,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param subx1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("subX1 to apply: %d\n", trolConfig->subx1);	


	/*************************** subx2 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->subx2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param subx1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("subX2 to apply: %d\n", trolConfig->subx2);

	/*************************** suby1 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->suby1,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param suby1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("suby1 to apply: %d\n", trolConfig->suby1);

	/*************************** suby2 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->suby2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param suby2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("suby2 to apply: %d\n", trolConfig->suby2);


	// CHECK SIZE IMAGE PARAMS

	if(trolConfig->subx2 > trolConfig->nx || trolConfig->subx1>trolConfig->subx2 || trolConfig->suby2 > trolConfig->ny || trolConfig->suby1 > trolConfig->suby2){
		printf("\n ERROR IN THE DIMENSIONS, PLEASE CHECK GIVEN VALUES \n ");
		exit(1);
	}
	if(trolConfig->subx1 == 0 && trolConfig->subx2==0 && trolConfig->suby1==0 && trolConfig->suby2==0){ // invert entire image 
		trolConfig->subx1 = 1;
		trolConfig->suby1 = 1;
		trolConfig->subx2 = trolConfig->nx;
		trolConfig->suby2 = trolConfig->ny;
	}

	return 1;

}



/**
 * 
 *  Method for read control file .
 *  @param fileParameters: name of file with control parameters. 
 *  @param trolConfig: structure to store configure information. 
 *  @param printlog: variable used to know if we must print log. 
 *  
 * */

int readTrolFile(char * fileParameters,  ConfigControl * trolConfig, int printLog){

	// try open the file with the 
	FILE * fReadParameters;
	float aux;
	char LINE [4096], * returnLine;
	char comment[200], name[100];
	fReadParameters = fopen(fileParameters, "r");
	if (fReadParameters == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;
	}
	int rfscanf; 
	
	/***************************  NUMBER OF CYCLES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->NumberOfCycles,comment);
	
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NumberOfCycles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	/*if( trolConfig->NumberOfCycles < 0){
		printf("milos: Error in NumberOfCycles parameter. review it. Not accepted: %d\n", trolConfig->NumberOfCycles);
		return 0;
	}*/
	//if(printLog) printf("Number Of Cycles  read from control file: %i\n", trolConfig->NumberOfCycles);
	if(printLog) printf("%s", LINE);

	/***************************  OBSERVED PROFILES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->ObservedProfiles,comment);
	if(trolConfig->ObservedProfiles[0]=='!')
		trolConfig->ObservedProfiles[0] = '\0';	
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Observed Profiles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Observed profiles file read from control file: %s\n", trolConfig->ObservedProfiles);
	if(printLog) printf("%s", LINE);

	/***************************  STRAY LIGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->StrayLightFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Stray light file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Stray light file read from control file: %s\n", trolConfig->StrayLightFile);
	if(printLog) printf("%s", LINE);

	/***************************  PSF FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->PSFFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param PSF file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	// first check if the input string is a valid file. 
	if(access(trolConfig->PSFFile,F_OK) == -1) { // is not a file
		double fwhm = atof(trolConfig->PSFFile);
		if(fwhm!=0){ // if it's cero then it isn't a correct value for fwhm and we will not use psf 
			
			trolConfig->FWHM = fwhm;
			trolConfig->ConvolveWithPSF = 1;
			if(printLog) printf("FWHM read from control file: %lf \n", trolConfig->FWHM);
		}
		else{
			trolConfig->ConvolveWithPSF = 0;
			trolConfig->FWHM = -1; // to indicate that we will not use FWHM 
		}
	} 
	else{
		
		if(strcmp(file_ext(trolConfig->PSFFile),PSF_FILE)==0 ){
			trolConfig->ConvolveWithPSF = 1;
		}
		else{
			printf("\nERROR: The extension of PSF file is not '.psf' . REVIEW IT %s, PLEASE\n",trolConfig->PSFFile);
			exit(EXIT_FAILURE);
		}

	}
	//if(printLog) printf("PSF file read from control file: %s. Convolve with PSF? %d \n", trolConfig->PSFFile,trolConfig->ConvolveWithPSF);
	if(printLog) printf("%s", LINE);

	/*************************** WAVELENGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->WavelengthFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param wavelength file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//	 check first if it's a file, second, is the extension is .grid or .fits . If it's a grid store in grid file, else store in wavelenghtfile 
	if(access(trolConfig->WavelengthFile,F_OK) != -1) { // is not a file
		if(strcmp(file_ext(trolConfig->WavelengthFile),GRID_FILE)==0){ // IT'S A GRID FILE
			strcpy(trolConfig->MallaGrid,trolConfig->WavelengthFile);
			strcpy(trolConfig->WavelengthFile,"");
			trolConfig->useMallaGrid = 1;
		}
		else if(strcmp(file_ext(trolConfig->WavelengthFile),FITS_FILE)==0){ // IT'S A FITS FILE
			strcpy(trolConfig->MallaGrid,"");
			trolConfig->useMallaGrid = 0;
		}
		else{
			printf("\n ERROR: The file given for the wavelengths does not have the correct extension. This must be .fits or .grid. Please check it. ");
			exit(EXIT_FAILURE);
		}
	}
	//if(printLog) printf("Wavelength grid file read from control file: %s\n", trolConfig->WavelengthFile);
	if(printLog) printf("%s", LINE);

	/*************************** ATOMIC PARAMETER  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->AtomicParametersFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Atomic parameters file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Atomic parameters file to apply: %s\n", trolConfig->AtomicParametersFile);
	if(printLog) printf("%s", LINE);

	/*************************** ABUNDANCES FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->AbundancesFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Abundances file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Abundances file to apply: %s\n", trolConfig->AbundancesFile);
	if(printLog) printf("%s", LINE);

	/*************************** INITIAL GUESS MODEL   FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel,comment);
	if(trolConfig->InitialGuessModel[0]=='!')
		trolConfig->InitialGuessModel[0] = '\0';
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Initial guess model 1 to apply: %s\n", trolConfig->InitialGuessModel);
	if(printLog) printf("%s", LINE);

	/*************************** INITIAL GUESS MODEL  2  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel_2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Initial guess model 2 to apply: %s\n", trolConfig->InitialGuessModel_2);
	if(printLog) printf("%s", LINE);

	/*************************** WEIGHT FOT STOKES I ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->WeightForStokes[0] = aux;
	//if(printLog) printf("Weight for Stokes I to apply: %lf\n", trolConfig->WeightForStokes[0]);
	if(printLog) printf("%s", LINE);

	/*************************** WEIGHT FOT STOKES Q ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes Q. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->WeightForStokes[1] = aux;
	//if(printLog) printf("Weight for Stokes Q to apply: %lf\n", trolConfig->WeightForStokes[1]);
	if(printLog) printf("%s", LINE);

	/*************************** WEIGHT FOT STOKES U ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes U. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->WeightForStokes[2] = aux;
	//if(printLog) printf("Weight for Stokes U to apply: %lf\n", trolConfig->WeightForStokes[2]);
	if(printLog) printf("%s", LINE);

	/*************************** WEIGHT FOT STOKES V ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes V. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->WeightForStokes[3] = aux;
	//if(printLog) printf("Weight for Stokes V to apply: %lf\n", trolConfig->WeightForStokes[3]);
	if(printLog) printf("%s", LINE);


	/*************************** AUTOMATIC SELECTED OF NODES **************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->automaticSelectOfNodes,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, automatic selected of nodes. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Automatic Selected of Nodes to apply: %i\n", trolConfig->automaticSelectOfNodes);
	if(printLog) printf("%s", LINE);


	///////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// READ PARAMETER FROM FIRST SPECTRAL LINE ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	/*************************** S0_SOURCE_FUNCTION_CONSTANT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[7],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for S_0 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for S_0 1 to apply: %d\n", trolConfig->fix[7]);
	if(printLog) printf("%s", LINE);

	/*************************** S1_SOURCE_FUNCTION_GRADIENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[8],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for S_1 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for S_1 1  to apply: %d\n", trolConfig->fix[8]);
	if(printLog) printf("%s", LINE);

	/*************************** ETA0_LINE_TO_CONTINUUM_ABSORPTION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[0],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for eta0 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for eta0 1 to apply: %d\n", trolConfig->fix[0]);
	if(printLog) printf("%s", LINE);

	/*************************** B_MAGNETIC_FIELD_STRENGTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[1],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for magnetic field 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for magnetic field 1 to apply: %d\n", trolConfig->fix[1]);
	if(printLog) printf("%s", LINE);

	/*************************** VLOS_LINE_OF_SIGHT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[2],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Line Of Sight Velocity. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Line Of Sight Velocity to apply: %d\n", trolConfig->fix[2]);
	if(printLog) printf("%s", LINE);

	/*************************** GM_MAGNETIC_FIELD_INCLINATION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[5],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for gamma 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for gamma 1 to apply: %d\n", trolConfig->fix[5]);
	if(printLog) printf("%s", LINE);

	/*************************** AZ_MAGNETIC_FIELD_AZIMUTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[6],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for phi 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for phi 1 to apply: %d\n", trolConfig->fix[6]);
	if(printLog) printf("%s", LINE);

	/*************************** DOPP_DOOPLER_WIDTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[3],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for lambda_doppler 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for lambda_doppler 1 to apply: %d\n", trolConfig->fix[3]);
	if(printLog) printf("%s", LINE);

	/*************************** AA_DAMPING_PARAMETER ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[4],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for damping 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for damping 1 to apply: %d\n", trolConfig->fix[4]);
	if(printLog) printf("%s", LINE);


	/*************************** MAC_MACROTURBULENT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[9],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Invert macroturbulence 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Invert macroturbulence 1 to apply: %d\n", trolConfig->fix[9]);
	if(printLog) printf("%s", LINE);


	///////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// READ PARAMETER FROM SECOND SPECTRAL LINE ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	/*************************** S0_SOURCE_FUNCTION_CONSTANT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[7],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for S_0 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for S_0 2 to apply: %d\n", trolConfig->fix2[7]);
	if(printLog) printf("%s", LINE);

	/*************************** S1_SOURCE_FUNCTION_GRADIENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[8],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for S_1 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for S_1 2  to apply: %d\n", trolConfig->fix2[8]);
	if(printLog) printf("%s", LINE);

	/*************************** ETA0_LINE_TO_CONTINUUM_ABSORPTION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[0],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for eta0 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for eta0 2 to apply: %d\n", trolConfig->fix2[0]);
	if(printLog) printf("%s", LINE);

	/*************************** B_MAGNETIC_FIELD_STRENGTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[1],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for magnetic field 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for magnetic field 2 to apply: %d\n", trolConfig->fix2[1]);
	if(printLog) printf("%s", LINE);

	/*************************** VLOS_LINE_OF_SIGHT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[2],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for LOS velocity 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for LOS velocity 2 to apply: %d\n", trolConfig->fix2[2]);
	if(printLog) printf("%s", LINE);

	/*************************** GM_MAGNETIC_FIELD_INCLINATION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[5],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for gamma 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for gamma 2 to apply: %d\n", trolConfig->fix2[5]);
	if(printLog) printf("%s", LINE);

	/*************************** AZ_MAGNETIC_FIELD_AZIMUTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[6],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for phi 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for phi 2 to apply: %d\n", trolConfig->fix2[6]);
	if(printLog) printf("%s", LINE);

	/*************************** DOPP_DOOPLER_WIDTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[3],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for lambda_doppler 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for lambda_doppler 2 to apply: %d\n", trolConfig->fix2[3]);
	if(printLog) printf("%s", LINE);

	/*************************** AA_DAMPING_PARAMETER ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[4],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Nodes for damping 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Nodes for damping 2 to apply: %d\n", trolConfig->fix2[4]);
	if(printLog) printf("%s", LINE);

	/*************************** MAC_MACROTURBULENT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix2[9],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Invert macroturbulence 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Invert macroturbulence 2 to apply: %d\n", trolConfig->fix2[9]);
	if(printLog) printf("%s", LINE);

	/*************************** STRAY LIGHT FACTOR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[10],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Stray light Factor. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->fix2[10] = trolConfig->fix[10];
	//if(printLog) printf("Stray light Factor  to apply: %d\n", trolConfig->fix[10]);
	if(printLog) printf("%s", LINE);


	/*************************** INVERT FILLING FACTOR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->InvertFillingFactor,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Filling Factor. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Filling Factor  to apply: %d\n", trolConfig->InvertFillingFactor);
	if(printLog) printf("%s", LINE);


	/*************************** MU ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->mu,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param mu=cos (theta). Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("mu=cos (theta) to apply: %f\n", trolConfig->mu);
	if(printLog) printf("%s", LINE);

	/*************************** EstimatedSNForI ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->EstimatedSNForI,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Estimated S/N for I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog){
		printf("%s", LINE);
		/*printf("Estimated S/N for I  to apply: %i\n", trolConfig->EstimatedSNForI);
		printf("Estimated noise for I  to apply: %lf\n", 1.0/trolConfig->EstimatedSNForI);*/
	} 
	trolConfig->noise = 1.0/trolConfig->EstimatedSNForI;
	trolConfig->sigma[0] = trolConfig->noise*trolConfig->noise;
	trolConfig->sigma[1] = trolConfig->sigma[0];
	trolConfig->sigma[2] = trolConfig->sigma[0];
	trolConfig->sigma[3] = trolConfig->sigma[0];
	// PUT VALUES IN ARRAY OF SIGMA 


 	/*************************** CONTINIUM CONTRAST ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ContinuumContrast,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Continuum contrast. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//trolConfig->toplim = trolConfig->ContinuumContrast;
	//if(printLog) printf("Continuum contrast to apply: %i\n", trolConfig->ContinuumContrast);
	if(printLog) printf("%s", LINE);


	/*************************** INITIAL_DIAGONAL_ELEMENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%f%99[^!]!",name, &aux,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial diagonal element. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	trolConfig->InitialDiagonalElement = aux;
	//if(printLog) printf("Initial diagonal element  to apply: %le\n", trolConfig->InitialDiagonalElement);
	if(printLog) printf("%s", LINE);
	
	/*************************** USE FFT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->useFFT,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Synthesis Profile Adjusted. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Save Synthesis Profile to apply: %d\n", trolConfig->SaveSynthesisAdjusted);
	if(printLog) printf("%s", LINE);

	return 1;

}


/**
 * 
 * Read init file for use with 
 * 
 * 
 * */
int readInitFile(char * fileParameters,  ConfigControl * trolConfig, int printLog){

	// try open the file with the 
	FILE * fReadParameters;
	char LINE [4096], * returnLine;
	char comment[200], name[100];
	fReadParameters = fopen(fileParameters, "r");
	if (fReadParameters == NULL)
	{
		if(printLog){
			printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
			printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		}
		return 0;
	}
	int rfscanf; 


	/***************************  TROL FILE  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->controlFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Control File. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Control File to apply: %s\n", trolConfig->controlFile);
	readTrolFile(trolConfig->controlFile,trolConfig,printLog);

	/***************************  type input stoke files ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->typeInputStokes,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Type Input Stokes. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Type Input Stokes to apply: %s\n", trolConfig->typeInputStokes);
	if(printLog) printf("%s", LINE);

	/***************************  type input stray light file  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->typeInputStrayLight,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param type input stray light. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Type input stray light to apply: %s\n", trolConfig->typeInputStrayLight);
	if(printLog) printf("%s", LINE);
	/*************************** NX ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->nx,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NX. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("NX to apply: %d\n", trolConfig->nx);
	if(printLog) printf("%s", LINE);
	/*************************** NY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ny,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NY. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("NY to apply: %d\n", trolConfig->ny);
	if(printLog) printf("%s", LINE);
	/*************************** subx1 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->subx1,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param subx1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("subX1 to apply: %d\n", trolConfig->subx1);	
	if(printLog) printf("%s", LINE);

	/*************************** subx2 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->subx2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param subx1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("subX2 to apply: %d\n", trolConfig->subx2);
	if(printLog) printf("%s", LINE);
	/*************************** suby1 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->suby1,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param suby1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("suby1 to apply: %d\n", trolConfig->suby1);
	if(printLog) printf("%s", LINE);
	/*************************** suby2 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->suby2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param suby2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("suby2 to apply: %d\n", trolConfig->suby2);
	if(printLog) printf("%s", LINE);
	// CHECK SIZE IMAGE PARAMS

	if(trolConfig->subx2 > trolConfig->nx || trolConfig->subx1>trolConfig->subx2 || trolConfig->suby2 > trolConfig->ny || trolConfig->suby1 > trolConfig->suby2){
		printf("\n ERROR IN THE DIMENSIONS, PLEASE CHECK GIVEN VALUES \n ");
		exit(1);
	}
	if(trolConfig->subx1 == 0 && trolConfig->subx2==0 && trolConfig->suby1==0 && trolConfig->suby2==0){ // invert entire image 
		trolConfig->subx1 = 1;
		trolConfig->suby1 = 1;
		trolConfig->subx2 = trolConfig->nx;
		trolConfig->suby2 = trolConfig->ny;
	}

	/***************************  output file prefix  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->outputPrefix,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param type input output prefix. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Output prefix to apply: %s\n", trolConfig->outputPrefix);
	if(printLog) printf("%s", LINE);

	/***************************  mask file prefix  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->MaskFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param maskfile. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Maskfile to apply: %s\n", trolConfig->MaskFile);
	if(printLog) printf("%s", LINE);
	/*************************** t1 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->t1,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param t1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("t1 to apply: %d\n", trolConfig->t1);
	if(printLog) printf("%s", LINE);
	/*************************** t2 ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->t2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param t2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("t2 to apply: %d\n", trolConfig->t2);
	if(printLog) printf("%s", LINE);
	/*************************** SAVE SYNTHESIS PROFILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->SaveSynthesisAdjusted,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Synthesis Profile Adjusted. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	//if(printLog) printf("Save Synthesis Profile to apply: %d\n", trolConfig->SaveSynthesisAdjusted);
	if(printLog) printf("%s", LINE);


	return 1;
}






/* Returns a pointer to the extension of 'string'.
 * If no extension is found, returns a pointer to the end of 'string'. */
char* file_ext(const char *string)
{
    assert(string != NULL);
    char *ext = strrchr(string, '.');
 
    if (ext == NULL)
        return (char*) string + strlen(string);
	 char * iter;	
    for (iter = ext + 1; *iter != '\0'; iter++) {
        if (!isalnum((unsigned char)*iter))
            return (char*) string + strlen(string);
    }
 
    return ext;
}

/**
 * 
 * 
 * 
 * */
char * get_basefilename (const char * fname) // returns the filename minus the extension
{
	char * ext;
	int i, j;
	ext = (char *)malloc(sizeof(char) * 4096);

	for ( i = strlen(fname)+1; i > 0; i--){
		if (fname[i] == '.'){
			for (j = 0; j < i; j++){
				ext[j] = fname[j];
			}
			ext[i] = '\0';
			i = 0;
		}
	}
	return ext;
}
