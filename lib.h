#include "defines.h"

/**
 * 
 */
int covarm(PRECISION *w, PRECISION *sig, float *spectro, int nspectro, PRECISION *spectra, PRECISION *d_spectra,
			  PRECISION *beta, PRECISION *alpha);




/**
 * 
 */
PRECISION *totalParcial(PRECISION * A, int f,int c,int dire);
/**
 * 
 */
PRECISION *totalParcialMatrix(PRECISION * A, int f,int c,int p);
/**
 * 
 */
PRECISION total(PRECISION * A, int f,int c);
/**
 * 
 */
int multmatrix(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col);
/**
 * 
 * */
int multmatrixCblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col);
/**
 * 
 */
PRECISION fchisqr(PRECISION * spectra,int nspectro,float *spectro,PRECISION *w,PRECISION *sig,PRECISION nfree);


/**
 * 
 */
int multmatrixIDLValue(PRECISION *a,int naf,int nac,PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);
/**
 * 
 */
void totalParcialMatrixf(PRECISION * A, int f,int c,int p,PRECISION *result);
/**
 * 
 */
void totalParcialf(PRECISION * A, int f,int c,PRECISION * result);
/**
 * 
 */
int multmatrix_transpose(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);

int multmatrix_transpose_param(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col);

int multmatrix_transpose_omp(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);

int multmatrix_transpose_cblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);


/**
 * 
 */
int CalculaNfree(int nspectro);


/**
 * 
 * */
PRECISION mean(PRECISION *dat, int numl);

/**
 * PRINT 
 * */
void printProgress (PRECISION percentage);

/*
* Check if path is a directory or not. 
*/
int isDirectory(const char *path);


/**
 * Own implementation of memcpy
 * */

void myMemCpy(PRECISION *dest, PRECISION *src, size_t n) ;



void strip_ext(char *fname);