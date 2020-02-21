
#include "defines.h"



/**
 * 
 * Clean the memory of the array "image"
 * @param image --> Array to clean memory 
 * @param numPixels --> Number of elements of array of image to clean. 
 */
void freeVpixels(vpixels * image, int numPixels);

/**
 * This function read the spectro image from the file "fitsFileSpectra" and store it into a struct of FitsImage
 * @param fitsFileSpectra --> name of the fits file to read 
 * Return the image read or NULL if something was wrong during the lecture. 
 */
FitsImage *  readFitsSpectroImage (const char * fitsFileSpectra, int forParallel);



/**
 * This function read the spectro image from the file "fitsFileSpectra" and store it into a struct of FitsImage
 * @param fitsFileSpectra --> name of the fits file to read 
 * Return the image read or NULL if something was wrong during the lecture. 
 */
FitsImage * readFitsSpectroNPixels (const char * fitsFileSpectra, int rowInit, int colInit,int rowEnd, int colEnd, int numPixelsRead,  int forParallel);

/**
 * This function read the lambda values for the image from the file "fitsFileLambda" and store it into a struct of FitsImage. The file of spectro must
 * be read it before call this method. 
 * @param fitsFileLambda --> name of the fits file to read with lambda values 
 * @param fitsImage --> struct of image 
 * Return 1 If the image has been read corectly if not return 0 
 */
int readFitsLambdaFile (const char * fitsFileLambda, FitsImage * fitsImage);

/**
 * 
 * */
PRECISION * readFitsLambdaToArray (const char * fitsFileLambda,  int * indexLine, int * nLambda);

/**
 * This function read the Stray Light values from the file "fitsFileStrayLight" and store it into a vector , in dimStrayLight will be stored the tam of fits file:
 *  (n lambdas or n lambdas X numPixels)
 * be read it before call this method. 
 * @param fitsFileLambda --> name of the fits file to read with lambda values 
 * @param fitsImage --> struct of image 
 * Return 1 If the image has been read corectly if not return 0 
 */
PRECISION * readFitsStrayLightFile (const char * fitsFileStrayLight, int * dimStrayLight, int numLambda);

/**
 * Clean the memory reserved for the image
 * @param image --> Image to clean memory . 
 */
void freeFitsImage(FitsImage * image);

/**
 * Write the models resutls in a fils file with 3 dimensiones: number of models x number of cols x number of rows. 
 * 
 * @param fitsFile --> Name of the file to store the image. 
 * @param numRows --> Number of rows of original image.
 * @param numCols --> Number of cols of original image.
 * @param vInitModel --> Array with the models obtained from the inversion. Each element of the array is a stucture with the models for one
 * pixel in the image. 
 * @param vChisqrf --> Array with the chisqrf calculated for each pixel in the image. 
 * @param fixed --> array with positions to write in the file, Positions are in the following order: 
 * [Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
 * @param addChisqr -->  to know if add chisqr at output model file 
 * 
 */
int writeFitsImageModels(const char * fitsFile, int numRows, int numCols, Init_Model * vInitModel, float * vChisqrf, int * vNumIterPixel, int addChiqr);
int writeFitsImageModelsSubSet(const char * fitsFileName, int numRows, int numCols,int rowInit, int colInit, int numPixelsWrite,int displs, Init_Model * vInitModel, float * vChisqrf, int * vNumIterPixel, int addChiqr);

int writeFitsImageModelsWithArray(char * fitsFile, int numRows, int numCols, PRECISION * eta0, PRECISION * B, PRECISION * vlos, PRECISION * dopp, PRECISION * aa, PRECISION * gm, PRECISION * az, PRECISION * S0, PRECISION * S1, PRECISION * mac, PRECISION * alfa, PRECISION * vChisqrf);
/**
 * Print the error status of status. 
 * @param status -> int code with the status error. 
 */
void printerror( int status);

/**
 * Write the image of profiles in a fits file with the same header of the fits spectro file. 
 * 
 * @param fitsProfileFile --> Name of the file to store the fits profile Image. 
 * @param fitsFileOrigin --> Name of the file with the origin spectro image. This image is used to copy the same header to the profile file. 
 * @param image --> Struct with the image procesed to store in the profile file. 
 * 
 */
int writeFitsImageProfiles(const char * fitsProfileFile, const char * fitsFileOrigin, FitsImage * image);


/**
 * 
 * This method is used to process the file with the parameters of the file to call the program MILOS. 
 * @param fileParameters 
 * @param maxIter 
 * @param clasicalEstimate
 * @param printSintesis
 * @param nameInputFileSpectra
 * @param nameOutputFileModels
 * @param nameOutputFilePerfiles
 * @param useConvolution
 * @param FWHM
 * @param DELTA
 * @param NMUESTRAS_G
 * 
 * @return 0 is there is something wrong reading the parameters, 1 is everything was fine. 
 */
/*int readParametersFileInput(char * fileParameters, int * maxIter, int * clasicalEstimate, int * printSintesis, char * nameInputFileSpectra, char * nameInputFileLambda,char * nameInputFileLines, char * nameInputFileInitModel,  PRECISION *  centralLambda, char *nameOutputFileModels, char * nameOutputFileProfiles, int * useConvolution, char * nameInputFilePSF, PRECISION * FWHM, int * KIND_CONVOLUTION);



int readFileCuanticLines(char * inputLineFile, PRECISION * cuanticDat, PRECISION centralLambda);

int readInitialModel(Init_Model * INIT_MODEL, char * fileInitModel);

int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, char * nameInputPSF);*/

/**
 * 
 * @param: fitsFile
 * @param: numRows
 * @param: numCols
 *  
 * */
int readSizeImageSpectro(const char * fitsFile, int * numRows, int * numCols);
