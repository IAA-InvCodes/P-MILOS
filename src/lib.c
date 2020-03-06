#include "defines.h"
#include <string.h>
#include "lib.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <sys/types.h>
#include <sys/stat.h>

//#include <omp.h>
//#include "mkl_cblas.h"



extern REAL AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
extern REAL * opa;


/*

 el tamaño de w es 	nlambda*NPARMS;

return 
	- beta de tam 1 x NTERMS
	- alpha de tam NTERMS x NTERMS

*/

int covarm(REAL *w,REAL *sig,float *spectro,int nspectro,REAL *spectra,REAL  *d_spectra,REAL *beta,REAL *alpha){	
	
	int j,i,bt_nf,bt_nc,aux_nf,aux_nc;

	
	REAL *BTaux,*APaux;

	//printf("\nVALORES DEL SIGMA SQUARE\n");

	for(j=0;j<NPARMS;j++){
		for(i=0;i<nspectro;i++){
			opa[i]=w[j]*(spectra[i+nspectro*j]-spectro[i+nspectro*j]);
		}

		BTaux=BT+(j*NTERMS);
		APaux=AP+(j*NTERMS*NTERMS);
		
		multmatrixIDLValue(opa,nspectro,1,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,BTaux,&bt_nf,&bt_nc,sig[j]); //bt de tam NTERMS x 1
		multmatrix_transpose(d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,APaux,&aux_nf,&aux_nc,w[j]/sig[j]);//ap de tam NTERMS x NTERMS
	}

	totalParcialf(BT,NPARMS,NTERMS,beta); //beta de tam 1 x NTERMS
	totalParcialMatrixf(AP,NTERMS,NTERMS,NPARMS,alpha); //alpha de tam NTERMS x NTERMS
	
	return 1;
}



REAL fchisqr(REAL * spectra,int nspectro,float *spectro,REAL *w,REAL *sig,REAL nfree){
	
	REAL TOT,dif;	
	REAL opa;
	int i,j;

	TOT=0;
	for(j=0;j<NPARMS;j++){
		opa=0;
		for(i=0;i<nspectro;i++){
			dif=spectra[i+nspectro*j]-spectro[i+nspectro*j];
			opa+= (dif*dif);
		}
		TOT+=((w[j]*opa)/(sig[j]));
		//TOT+=((w[j]*opa)/(sig[j]*sig[j]));
		//TOT+= opa;///(sig[j]*sig[j]);
	}
		
	//return TOT/15;		
	return TOT/nfree;
	
}


/*

	Multiplica la matriz a (tamaño naf,nac)
	por la matriz b (de tamaño nbf,nbc)
	al estilo IDL, es decir, filas de a por columnas de b,
	el resultado se almacena en resultOut (de tamaño fil,col)

	El tamaño de salida (fil,col) corresponde con (nbf,nac).

	El tamaño de columnas de b, nbc, debe de ser igual al de filas de a, naf.

*/
int multmatrixIDLValue(REAL *a,int naf,int nac,REAL *b,int nbf,int nbc,REAL *result,int *fil,int *col,REAL value){
    
   int i,j,k;
   REAL sum;
	
	if(naf==nbc){
		(*fil)=nbf;
		(*col)=nac;

		for ( i = 0; i < nbf; i++){
		    for ( j = 0; j < nac; j++){
				sum=0;
				for ( k = 0;  k < naf; k++){
					//printf("i: %d,j:%d,k=%d .. a[%d][%d]:%f  .. b[%d][%d]:%f\n",i,j,k,k,j,a[k*nac+j],i,k,b[i*nbc+k]);
					sum += a[k*nac+j] * b[i*nbc+k];
				}
				//printf("Sum, result[%d][%d] : %f \n",i,j,sum);
				result[((nac)*i)+j] = sum/value;
      		} 
		}
		return 1;
	}
	else
		printf("\n \n Error en multmatrixIDLValue no coinciden nac y nbf!!!! ..\n\n");
	return 0;
}


void totalParcialf(REAL * A, int f,int c,REAL * result){

	int i,j;
	REAL sum;
	for(i=0;i<c;i++){
		sum = 0;
		for(j=0;j<f;j++){
			sum+=A[j*c+i];
		}
		result[i] = sum;
	}
}


void totalParcialMatrixf(REAL * A, int f,int c,int p,REAL *result){

	int i,j,k;
	REAL sum;
	for(i=0;i<f;i++)
		for(j=0;j<c;j++){
			sum=0;
			for(k=0;k<p;k++)
				sum+=A[i*c+j+f*c*k];
			result[i*c+j] = sum;
		}

//	return result;
}

/*
	Multiplica la matriz a (tamaño naf,nac)
	por la matriz b (de tamaño nbf,nbc)
	al estilo multiplicación algebraica de matrices, es decir, columnas de a por filas de b,
	el resultado se almacena en resultOut (de tamaño fil,col)

	El tamaño de salida (fil,col) corresponde con (nbf,nac).

	El tamaño de columnas de a, nac, debe de ser igual al de filas de b, nbf.
*/

int multmatrix(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col){
    
    int i,j,k;
    PRECISION sum;
    
	if(nac==nbf){
		(*fil)=naf;
		(*col)=nbc;
		
//		free(*result);
//		*result=calloc((*fil)*(*col),sizeof(PRECISION));
//		printf("a ver ..\n");

		for ( i = 0; i < naf; i++)
		    for ( j = 0; j < nbc; j++){
				sum=0;
				for ( k = 0;  k < nbf; k++){
//					printf("i: %d,j:%d,k=%d .. a[%d][%d]  .. b[%d][%d]\n",i,j,k,i,k,k,j);
					sum += a[i*nac+k] * b[k*nbc+j];
				}
//				printf("Sum\n");
				result[(*col)*i+j] = sum;

      		} 

		return 1;
	}
	return 0;

}




int multmatrix_transpose(REAL *a,int naf,int nac, REAL *b,int nbf,int nbc,REAL *result,int *fil,int *col,REAL value){
    
    int i,j,k;
    REAL sum;
    
	if(nac==nbc){
		(*fil)=naf;
		(*col)=nbf;
		
		for ( i = 0; i < naf; i++){
		    for ( j = 0; j < nbf; j++){
				sum=0;
				for ( k = 0;  k < nbc; k++){
					sum += a[i*nac+k] * b[j*nbc+k];
				}

				result[(*col)*i+j] = (sum)*value;
     		} 
		
		}
		return 1;
	}else{
		printf("\n \n Error en multmatrix_transpose no coinciden nac y nbc!!!! ..\n\n");
	}

	return 0;
}


//Media de un vector de longitud numl
PRECISION mean(PRECISION *dat, int numl){
	
	PRECISION auxsum;
	int i;

	auxsum=0;
	for(i=0;i<numl;i++){
		auxsum=auxsum+dat[i];		
	}

	return auxsum/numl;
}

/**
 * 
 */
int CalculaNfree(int nspectro)
{
	int nfree;
	nfree = 0;

	nfree = (nspectro * NPARMS) - NTERMS;

	return nfree;
}



void printProgress (PRECISION percentage)
{
    //int val = (int) (percentage * 100);
    //int lpad = (int) (percentage * PBWIDTH);
    //int rpad = PBWIDTH - lpad;
    //printf ("\r%3f%% [%.*s%*s]",  percentage, lpad, PBSTR, rpad, "");
	 printf ("\r%3f %%",  percentage);
    fflush (stdout);
}


int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}


void myMemCpy(PRECISION *dest, PRECISION *src, size_t n) 
{ 
   // Typecast src and dest addresses to (char *) 
	size_t i;
   // Copy contents of src[] to dest[] 
   for (i=0; i<n; i++) 
      dest[i] = src[i]; 
} 



void strip_ext(char *fname)
{
    char *end = fname + strlen(fname);

    while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }
    if ((end > fname && *end == '.') &&
        (*(end - 1) != '\\' && *(end - 1) != '/')) {
        *end = '\0';
    }  
}



