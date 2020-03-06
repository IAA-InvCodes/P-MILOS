
#include "defines.h"

//int readConfigControl(char * configFile, ConfigControl * trolConfig, int printLog);

PRECISION readFileCuanticLines(const char * inputLineFile, PRECISION * cuanticDat, int line2Read, int printLog);

//int readInitialModel(Init_Model * INIT_MODEL, const char * fileInitModel);

int readInitialModel(Init_Model * INIT_MODEL, char * fileInitModel);

int readMallaGrid(const char * fileMallaGrid, PRECISION * initialLambda, PRECISION * step, PRECISION * finalLambda, int printLog);

int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, const char * nameInputPSF, PRECISION centralWaveLenght);

void loadInitialValues(ConfigControl * configControlFile);

int readParametersFileInput(char * fileParameters,  ConfigControl * trolConfig, int printLog);

int readTrolFile(char * fileParameters,  ConfigControl * trolConfig, int printLog);

int readInitFile(char * fileParameters,  ConfigControl * trolConfig, int printLog);

char* file_ext(const char *string);

char * get_basefilename (const char * fname);
