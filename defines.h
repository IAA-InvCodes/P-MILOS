
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "fitsio.h"
//#include <fftw3.h> //siempre detras de complex.h!
//#include <math.h>
//#include <stdio.h>



#ifndef DEFINES_H_
#define DEFINES_H_

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
// USER CONFIGURATION

//#define CENTRAL_WL   6173.341000 //6173.341000 (oficial requ.) //6173.3500// 341000  // // 6173.335600 //6173.335400


//NumeroS cuanticos
#define CUANTIC_NWL 1
#define CUANTIC_SLOI 2
#define CUANTIC_LLOI 1
#define CUANTIC_JLOI 1
#define CUANTIC_SUPI 2
#define CUANTIC_LUPI 2
#define CUANTIC_JUPI 0


#define NOISE_SIGMA 0.001 

#define CLASSICAL_ESTIMATES_SAMPLE_REF 4 //Muestra referencia para cambio de cuadrante de azimuth. Depende del numero de muestras y posicion Continuo

#define NTERMS 11  //ojo si es mayor q 10 casca el svdCordic (esta version)

//##############################################
//SVD CONFIGURATION
#define USE_SVDCMP 1    //1 for using SVDCMP and 0 for using SVD_CORDIC  -->  Note: the SVDCMP doesn't work in float! only double

#define NORMALIZATION_SVD 0 //1 for using normalization matrixes ONLY  in the SVD_CORDIC

#define NUM_ITER_SVD_CORDIC 36 //9,18,27,36  --> 18 parece ok!

#define LIMITE_INFERIOR_PRECISION_SVD pow(2.0,-54)
#define LIMITE_INFERIOR_PRECISION_TRIG pow(2.0,-54)
#define LIMITE_INFERIOR_PRECISION_SINCOS pow(2.0,-54)
#define PRECISION double //double or float

//#############################################



// END USER CONFIGURATION
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------

// DONT'T MODIFY ANYTHING BELOW OF THIS LINE

#define PI 	3.14159265358979323846264338327950288419716939937510
 		 	
#define ILAMBDA 0.1
#define TOPLIM 1e-12 
#define SLIGHT 0
#define NOISE 1e-10 //0.001

#define RR  0.5641895836

#define VLIGHT 2.99792458e+5 //;light speed (km/s); 

#define CTE4_6_13 4.6686411e-13
#define AH 1.0 //angulo heliocentrico

#define FFT_FORWARD -1 
#define FFT_BACKWARD +1

#define NPARMS 4 //(IQUV)

#define LONG_PUNTERO_CALCULOS_COMPARTIDOS 100  //no se usa...

#define INSTRUMENTAL_CONVOLUTION_INTERPOLACION 0  //realizar interpolacion en la convolucion ?? //No funciona !

//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
struct INIT_MODEL{
	PRECISION eta0; // 0
	PRECISION B;//magnetic field    
	PRECISION vlos;
	PRECISION dopp;
	PRECISION aa;
	PRECISION gm; //5
	PRECISION az;
	PRECISION S0;
	PRECISION S1;
	PRECISION mac; //9
	PRECISION alfa;		
};

struct CUANTIC{  
	
	PRECISION N_PI;
	PRECISION N_SIG;
	PRECISION * NUB;//size stored in  n_sig
	PRECISION * NUP;//size stored in n_pi
	PRECISION * NUR;//size stored in n_sig
	PRECISION * WEB;//size stored in n_sig
	PRECISION * WEP;//size stored in n_pi
	PRECISION * WER;//size stored in n_sig
	PRECISION GL;
	PRECISION GU;
	PRECISION GEFF;
	PRECISION FO;	
	
};

typedef struct INIT_MODEL Init_Model;
typedef struct CUANTIC Cuantic;

/******************************************************/

void InitializePointerShareCalculation();
void ResetPointerShareCalculation();
void ReadPointerShareCalculation(int Numero,PRECISION ** a,...);
void AsignPointerShareCalculation(int Numero,PRECISION * a,...);
void DeleteSpectraCalculation();

void AllocateMemoryDerivedSynthesis(int numl);
void FreeMemoryDerivedSynthesis();

/******************************************************/


Cuantic * create_cuantic(PRECISION * dat, int log);

int me_der(Cuantic *cuantic,Init_Model *initModel,PRECISION * wlines,PRECISION *lambda,int nlambda,
			PRECISION *d_spectraOut,PRECISION *spectra, PRECISION * spectra_slight,PRECISION ah,PRECISION * slight,int calcSpectra,int filter);

int mil_sinrf(Cuantic *cuantic,Init_Model *initModel,PRECISION * wlines,PRECISION *lambda,int nlambda,PRECISION *spectra,
			PRECISION ah,PRECISION * slight,PRECISION * spectra_mc, int filter);
			

PRECISION * fgauss(PRECISION MC, PRECISION * eje,int neje,PRECISION landa,int deriv);
PRECISION * fgauss_WL(PRECISION FWHM, PRECISION step_between_lw, PRECISION lambda0, PRECISION lambdaCentral, int nLambda, int * sizeG);



int Guarda(char * nombre,PRECISION *v,int nv);
int GuardaC(char * nombre,PRECISION _Complex *v,int nv,int a);

int fvoigt(PRECISION damp,PRECISION *vv,int nvv,PRECISION *h, PRECISION *f);

//PRECISION * vgauss(PRECISION fwhm,int nmuestras_G,PRECISION delta);

/******************* DEFINITIONS FOR READ FITS FILE *********************/

/* 
	Every fits_image will be store in memory like  4 Vector of PRECISION TYPE using this dimension: rows*cols*nLambdas. 
	If we want access a pixel in the vector: stockesI[ (numRow*numCol) + (rows*cols*nLamba)] (been lambda from 0 to nLambda-1)
*/



#define CTYPE1 "CTYPE1"
#define CTYPE2 "CTYPE2"
#define CTYPE3 "CTYPE3"
#define CTYPE4 "CTYPE4"
#define CUNIT1 "CUNIT1"
#define CUNIT2 "CUNIT2"
#define CUNIT3 "CUNIT3"

#define CUNIT_ANSTROM "Angstrom"
#define CUNIT_ARCSEC "arcsec"
#define CTYPE_WAVE "WAVE-GRI"
#define CTYPE_HPLN_TAN "HPLN-TAN"
#define CTYPE_HPLT_TAN "HPLT-TAN"
#define CTYPE_STOKES "STOKES"


/* This Sttructure will store information relative a pixel in the image to process by fuction lmils in C, MPI and CUDA. 
	In order to do more flexible the structure and pass this structure in the message to MPI and KERNEL of CUDA, we have decided 
	put an attribute for number of pixels to store in spectro. We will play with this parameter to pass more or less pixels through MPI and to the 
	KERNEL of CUDA. 
 */
struct VPIXEL {
	PRECISION * vLambda;
	float * spectro;
	int nLambda;
};

typedef struct VPIXEL vpixels;



struct FITS_IMAGE{  

	int rows;  // number of rows in the image
	int cols;  // number of cols in the image
	int nLambdas; // number of lambdas in the image 
	int numStokes; // number of stokes paramters, normally is 4

	/* POSITION OF EACH DATA IN THE DIMENSIONS OF FITS */
	int pos_lambda;
	int pos_row;
	int pos_col;
	int pos_stokes_parameters; 

	char ctype_1[FLEN_CARD];
	char ctype_2[FLEN_CARD];
	char ctype_3[FLEN_CARD];
	char ctype_4[FLEN_CARD];
	char cunit_1[FLEN_CARD];
	char cunit_2[FLEN_CARD];
	char cunit_3[FLEN_CARD];
	char cunit_4[FLEN_CARD];

	int numPixels;
	vpixels * pixels;
	PRECISION * vLambdaImagen;
	float * spectroImagen;
};

typedef struct FITS_IMAGE FitsImage;

#define NUMBER_PARAM_MODELS 11

#define PATH_MAX 4096

struct CONFIG_CONTROL{

	int NumberOfCycles;
	char ObservedProfiles[4096];
	char StrayLightFile[4096];
	char PSFFile[4096];
	char WavelengthFile[4096];
	char AtomicParametersFile[4096];
	char InitialGuessModel[4096];
	char InitialGuessModel_2[4096];
	PRECISION WeightForStokes[4];
	int InvertFillingFactor;
	int InvertStrayLightFactor;
	PRECISION mu;
	int EstimatedSNForI;
	int ContinuumContrast;
	PRECISION ToleranceForSVD;
	PRECISION InitialDiagonalElement;
	int ConvolveWithPSF;
	PRECISION FWHM;
	PRECISION CentralWaveLenght;
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	int fix[11]; // eta0, B , vlos, dopp, aa, gm , az, S0, S1, mac, alpha
	int fix2[11]; // eta0, B, vlos, dopp, aa, gm , az, S0, S1, mac, alpha
	int saveChisqr;
	PRECISION toplim; // Optional minimum relative difference between two succesive merit-function values
	PRECISION sigma [4];
	PRECISION noise;
	int UseClassicalEstimates;
	int UseRTEInversion;
	int SaveSynthesisAdjusted;	
	int typeFileOutputModel; // 0 print to FITS , 1 print to TXT. 
	char OutputModelFile[4096];
	char OutputSynthesisFile[4096];
	char MallaGrid[4096];
	char AbundancesFile[4096];
	int useMallaGrid; // value 1 --> use malla grid, value 0 --> use fits fil
	int automaticSelectOfNodes;
	char controlFile [4096];
	char typeInputStokes [50];
	char typeInputStrayLight[50];
	int nx;
	int ny;
	int subx1;
	int subx2;
	int suby1;
	int suby2;	
	char outputPrefix[4096];	
	char MaskFile[4096];
	int t1;
	int t2;
};

typedef struct CONFIG_CONTROL ConfigControl;

// CONSTANTS TO COMPARE DATA FOR READ FROM CONFIGURATION FILE 
#define NUMBER_OF_CYCLES "NumberOfCycles"
#define OBSERVED_PROFILES "ObservedProfiles"
#define STRAY_LIGHT_FILE "StrayLightFile"
#define PSF_FILE "PSFFile"
#define WAVE_LENGHT_FILE "WavelengthFile"
#define ATOMIC_PARAMETERS_FILE "AtomicParametersFile"
#define INITIAL_GUESS_MODEL "InitialGuessModel"
#define WEIGHT_FOR_STOKESI  "WeightForStokesI"
#define WEIGHT_FOR_STOKESQ  "WeightForStokesQ"
#define WEIGHT_FOR_STOKESU  "WeightForStokesU"
#define WEIGHT_FOR_STOKESV  "WeightForStokesV"
#define INVERT_MACROTURBULENCE "InvertMacroturbulence"
#define INVERT_FILLING_FACTOR "InvertFillingFactor"
#define INVERT_STRAY_LIGHT_FACTOR "InvertStrayLightFactor"
#define MU "mu"
#define ESTIMATEDSNFORI "EstimatedSNForI"
#define CONTINUUM_CONTRAST "ContinuumContrast"
#define TOLERANCE_FOR_SVD "ToleranceForSVD"
#define INITIAL_DIAGONAL_ELEMENT "InitialDiagonalElement"
#define USE_INTERPOLAR_SPLINES_OR_LINEAR "UseInterpolarSplinesOrLinear"
#define CONVOLVE_WITH_PSF "ConvolveWithPSF"
#define FWHM_FILE "FWHM"
#define TYPE_CONVOLUTION "TypeConvolution"
#define GAS_PRESSURE_AT_SURFACE_1 "GasPressureAtSurface1" 
#define GAS_PRESSURE_AT_SURFACE_2 "GasPressureAtSurface2" 
#define MAGNETIC_PRESSURE_TERM "MagneticPressureTerm"
#define NTL "ntl"
#define NLIOBS "nliobs"
#define CENTRAL_WAVE_LENGHT "CentralWaveLenght"
#define ETA0_LINE_TO_CONTINUUM_ABSORPTION "ETA0_LineToContiniuumAbsorption"
#define B_MAGNETIC_FIELD_STRENGTH "B_MagneticFieldStrength"
#define VLOS_LINE_OF_SIGHT_VELOCITY "VLOS_LineOfSightVelocity"
#define DOPP_DOOPLER_WIDTH "DOPP_DooplerWidth"
#define AA_DAMPING_PARAMETER "AA_DampingParameter"
#define GM_MAGNETIC_FIELD_INCLINATION "GM_MagneticFieldInclination"
#define AZ_MAGNETIC_FIELD_AZIMUTH "AZ_MagneticFieldAzimuth"
#define S0_SOURCE_FUNCTION_CONSTANT "S0_SourceFunctionConstant"
#define S1_SOURCE_FUNCTION_GRADIENT "S1_SourceFunctionGradient"
#define MAC_MACROTURBULENT_VELOCITY "MAC_MacroturbulentVelocity"
#define ALPHA_FILLING_FACTOR "ALPHA_FillingFactor"
#define SAVE_CHISQR "SaveChisqr"
#define USE_CLASSICAL_ESTIMATES "UseClassicalEstimates"
#define USE_RTE_INVERSION "UseRTEInversion"
#define SAVE_SYNTHESIS_PROFILE "SaveSynthesisAdjusted"
#define OUTPUT_MODEL_FILE "OutputModelFile"
#define OUTPUT_SYNTHESIS_FILE "OutputSynthesisFile"
#define SIGMA_FILE "sigma"
#define NOISE_FILE "noise"
#define TOPLIM_FILE "toplim" 

// values for init guess model 

#define INITIAL_MODEL_ETHA0 "INITIAL_MODEL_ETHA0"
#define INITIAL_MODEL_B "INITIAL_MODEL_B"
#define INITIAL_MODEL_VLOS "INITIAL_MODEL_VLOS"
#define INITIAL_MODEL_LAMBDADOPP "INITIAL_MODEL_LAMBDADOPP"
#define INITIAL_MODEL_AA "INITIAL_MODEL_AA"
#define INITIAL_MODEL_GM "INITIAL_MODEL_GM"
#define INITIAL_MODEL_AZI "INITIAL_MODEL_AZI"
#define INITIAL_MODEL_S0 "INITIAL_MODEL_S0"
#define INITIAL_MODEL_S1 "INITIAL_MODEL_S1"
#define INITIAL_MODEL_MAC "INITIAL_MODEL_MAC"
#define INITIAL_MODEL_ALFA "INITIAL_MODEL_ALFA"


// TYPES OF CONVOLUTION 

#define CONVOLUTION_FFT "FFT"
#define CONVOLUTION_DIRECT "DIRECT"

// PERCENTAGE OF ITERATION
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

struct NAME_FILE {
	char  name [PATH_MAX];
};

typedef struct NAME_FILE nameFile;



#define PER_FILE ".per"
#define GRID_FILE ".grid"
#define FITS_FILE ".fits"
#define TROL_FILE ".mtrol"
#define MOD_FILE ".mod"
#define OUTPUT_MOD_FIT_EXT "_output.fits"
#define MOD_FITS "_mod.fits"
#define STOKES_FIT_EXT "_stokes.fits"
#define OUTPUT_MOD_TXT_EXT "_output_mod.txt"
#define STOKES_PER_EXT "_stokes.per"


#endif /*DEFINES_H_*/


