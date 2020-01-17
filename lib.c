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



extern PRECISION AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
extern PRECISION * opa;


/*

 el tamaño de w es 	nlambda*NPARMS;

return 
	- beta de tam 1 x NTERMS
	- alpha de tam NTERMS x NTERMS

*/

int covarm(PRECISION *w,PRECISION *sig,float *spectro,int nspectro,PRECISION *spectra,PRECISION  *d_spectra,
		PRECISION *beta,PRECISION *alpha){	
	
	int j,i,bt_nf,bt_nc,aux_nf,aux_nc;
	//static PRECISION AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
	//PRECISION opa[nspectro];
	
	PRECISION *BTaux,*APaux;
	PRECISION auxWeight,aux;
	//PRECISION *opa;
	//opa = calloc(nspectro,sizeof(PRECISION));
	//printf("\nVALORES DEL SIGMA SQUARE\n");

	for(j=0;j<NPARMS;j++){
		for(i=0;i<nspectro;i++){
			opa[i]=w[j]*(spectra[i+nspectro*j]-spectro[i+nspectro*j]);
		}
		auxWeight = (w[j]/(sig[j]));
		//auxWeight = (w[j]);

		BTaux=BT+(j*NTERMS);
		APaux=AP+(j*NTERMS*NTERMS);
		
		multmatrixIDLValue(opa,nspectro,1,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,BTaux,&bt_nf,&bt_nc,sig[j]); //bt de tam NTERMS x 1
		//multmatrix_transpose(d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,APaux,&aux_nf,&aux_nc,w[j]/(sig[j]*sig[j]));//ap de tam NTERMS x NTERMS
		multmatrix_transpose(d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,APaux,&aux_nf,&aux_nc,auxWeight);//ap de tam NTERMS x NTERMS
		//multmatrix_transpose_param(d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,d_spectra+j*nspectro*NTERMS,NTERMS,nspectro,APaux,&aux_nf,&aux_nc);//ap de tam NTERMS x NTERMS
		//printf("\n FIRST PART OF ALFA MATRIX:\n ");
		//for(i=0;i<aux_nf*aux_nc;i++){
			//printf("\nVALUE OF APAUX[I]=%lf  \t VALUE OF AUX WEIGHT=%lf",APaux[i], auxWeight);
			//APaux[i] *= auxWeight;
			/*aux = APaux[i] * auxWeight;
			APaux[i] =  aux;*/
			//printf("\n VALUE OF APAUX[i] alter mult: %lf",APaux[i]);
		//}
	}

	totalParcialf(BT,NPARMS,NTERMS,beta); //beta de tam 1 x NTERMS
	totalParcialMatrixf(AP,NTERMS,NTERMS,NPARMS,alpha); //alpha de tam NTERMS x NTERMS
	
	return 1;
}



PRECISION fchisqr(PRECISION * spectra,int nspectro,float *spectro,PRECISION *w,PRECISION *sig,PRECISION nfree){
	
	PRECISION TOT,dif;	
	PRECISION opa;
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
int multmatrixIDLValue(PRECISION *a,int naf,int nac,PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value){
    
     int i,j,k;
    PRECISION sum;
	
	if(naf==nbc){
		(*fil)=nbf;
		(*col)=nac;
		
//		free(*result);
//		result=calloc((nbf)*(nac),sizeof(PRECISION));
//		printf("a ver ..\n");

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


/*
	dire:
		1: suma por filas, return PRECISION * de tam f
		2: suma por columnas, return PRECISION * de tam c
*/

PRECISION *totalParcial(PRECISION * A, int f,int c,int dire){

	int i,j;
//	PRECISION 	sum;
	PRECISION *result;	
	result=calloc(dire==1?f:c,sizeof(PRECISION));

	for(i=0;i<f;i++)
		for(j=0;j<c;j++){
			result[(dire==1)?i:j]+=A[i*c+j];
		}

	return result;
}

void totalParcialf(PRECISION * A, int f,int c,PRECISION * result){

	int i,j;
//	PRECISION 	sum;

//	result=calloc(dire==1?f:c,sizeof(PRECISION));
	
	for(i=0;i<c;i++){
		result[i]=0;
		for(j=0;j<f;j++){
			result[i]+=A[j*c+i];
		}
	}
}


/*
return matriz de tam f*c
*/

PRECISION *totalParcialMatrix(PRECISION * A, int f,int c,int p){

	int i,j,k;
//	PRECISION 	sum;
	PRECISION *result;	
	result=calloc(f*c,sizeof(PRECISION));

	for(i=0;i<f;i++)
		for(j=0;j<c;j++){
			for(k=0;k<p;k++)
				result[i*c+j]+=A[i*c+j+f*c*k];
		}

	return result;
}

void totalParcialMatrixf(PRECISION * A, int f,int c,int p,PRECISION *result){

	int i,j,k;
//	PRECISION 	sum;
//	PRECISION *result;	
//	result=calloc(f*c,sizeof(PRECISION));

	for(i=0;i<f;i++)
		for(j=0;j<c;j++){
			result[i*c+j]=0;
			for(k=0;k<p;k++)
				result[i*c+j]+=A[i*c+j+f*c*k];
		}

//	return result;
}


PRECISION total(PRECISION * A, int f,int c){

	int i,j;
	PRECISION 	sum;
	sum=0;
	for(i=0;i<f;i++)
		for(j=0;j<c;j++)
			sum+=A[i*c+j];

	return sum;
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


/*
	Multiplica la matriz a (tamaño naf,nac)
	por la matriz b (de tamaño nbf,nbc)
	al estilo multiplicación algebraica de matrices, es decir, columnas de a por filas de b,
	el resultado se almacena en resultOut (de tamaño fil,col)

	El tamaño de salida (fil,col) corresponde con (nbf,nac).

	El tamaño de columnas de a, nac, debe de ser igual al de filas de b, nbf.
*/

int multmatrixCblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col){
    
    
    
    
	if(nac==nbf){
		(*fil)=naf;
		(*col)=nbc;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,naf, nbc, nac, 1.0, a, nac, b, nbc, 0.0, result, nbc);
		
		return 1;
	}
	return 0;

}



int multmatrix_transpose(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value){
    
    int i,j,k;
    PRECISION sum;
    
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

int multmatrix_transpose_param(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col){
    
    int i,j,k;
    PRECISION sum;
    
	if(nac==nbc){
		(*fil)=naf;
		(*col)=nbf;
		
		for ( i = 0; i < naf; i++){
		    for ( j = 0; j < nbf; j++){
				sum=0;
				for ( k = 0;  k < nbc; k++){
					sum += a[i*nac+k] * b[j*nbc+k];
				}

				result[(*col)*i+j] = (sum);
     		} 
		
		}
		return 1;
	}else{
		printf("\n \n Error en multmatrix_transpose no coinciden nac y nbc!!!! ..\n\n");
	}

	return 0;
}




/*int multmatrix_transpose_omp(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value){
    

	
	if(nac==nbc){
		(*fil)=naf;
		(*col)=nbf;
		//PRECISION matrix_aux[nac*naf];
		#pragma omp parallel
		{
			int i,j,k;
			//PRECISION sum;
			#pragma omp for 
			for ( i = 0; i < naf; i++){
				for ( j = 0; j < nbf; j++){
					result[(*col)*i+j]=0;
					for ( k = 0;  k < nbc; k++){
						result[(*col)*i+j] += a[i*nac+k] * b[j*nbc+k];
					}
					result[(*col)*i+j]=result[(*col)*i+j]*value;
				} 
			}


		}
		return 1;
	}else{
		printf("\n \n Error en multmatrix_transpose no coinciden nac y nbc!!!! ..\n\n");
	}

	return 0;
}
*/

/**
 * Mul matrix using cblas from GSL
 * 
 * */
int multmatrix_transpose_cblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value){


	
	if(nac==nbc){
		
		(*fil)=naf;
		(*col)=nbf;

		/*gsl_matrix_view matrixA = gsl_matrix_view_array(a, naf, nac);
   	gsl_matrix_view matrixB = gsl_matrix_view_array(b, nbf, nbc);
   	gsl_matrix_view matrixC = gsl_matrix_view_array(result, naf, nbf);
  		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, value , &matrixA.matrix, &matrixB.matrix, 0.0, &matrixC.matrix);*/
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,naf, nbc, nac, 1.0, a, nbc, b, nbc, value, result, nbc);
		// as is ROW MAJOR AND NO TRANSPOSE THE MULTIPLICATION, THEN LDA = NAC , LDB = NBC and LDC = nbc
		//#pragma omp parallel
		{
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,naf, nbf, nac, value, a, nac, a, nac, 0.0, result, naf);
		}
		
		/*gsl_matrix_view matrixResult = gsl_matrix_view_array(result,nac,nbf);
		*/

		return 1;
	}else{
		printf("\n \n Error en multmatrix_transpose no coinciden naf %d nac %d y nbf %d nbc %d !!!! .. \n\n", naf, nac, nbf, nbc);
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



/*
 * Cambiamos la forma para tener en cada fila I Q U V
 * Tambien reajusta el tamaño para eliminar las posiciones vacias
 */

void reformarVector(PRECISION **spectro,int neje){
	
	PRECISION *aux;
	int i;
	aux=(PRECISION *)calloc(neje*4,sizeof(PRECISION));

	for(i=0;i<neje;i++){
		aux[i]=(*spectro)[i*4];		
		aux[neje+i]=(*spectro)[i*4+1];
		aux[2*neje+i]=(*spectro)[i*4+2];		
		aux[3*neje+i]=(*spectro)[i*4+3];		
	}
	
	free(*spectro);
	*spectro=aux;
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