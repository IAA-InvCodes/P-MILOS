#include "utilsFits.h"
#include "fitsio.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <locale.h>
#include <unistd.h>

/**
 * Clean memory from fits image
 */
void freeVpixels(vpixels * image, int numPixels){
	int i;
	for(i=9;i<numPixels;i++){
		free(image[i].spectro);
		free(image[i].vLambda);
	}
	free(image);
}


/**
 * This function read the spectro image from the file "fitsFileSpectra" and store it into a struct of FitsImage
 * fitsFileSpectra --> name of the fits file to read 
 * Return the image read or NULL if something was wrong during the lecture. 
 */
/*FitsImage *  readFitsSpectroImage (const char * fitsFileSpectra, int forParallel){
   fitsfile *fptr;   
   FitsImage * image =  malloc(sizeof(FitsImage));
   int status = 0;   
	PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, naxis, anynul, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; 
	char comment[FLEN_CARD];   
	
	int i, j, k, h;
   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   if (!fits_open_file(&fptr, fitsFileSpectra, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;


			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_lambda; 
			int pos_row;
			int pos_col;
			int pos_stokes_parameters;
			// LAMBDA POSITION
			if(strcmp(image->ctype_1,CTYPE_WAVE)==0) pos_lambda = 0;
			if(strcmp(image->ctype_2,CTYPE_WAVE)==0) pos_lambda = 1;
			if(strcmp(image->ctype_3,CTYPE_WAVE)==0) pos_lambda = 2;
			if(strcmp(image->ctype_4,CTYPE_WAVE)==0) pos_lambda = 3;

			// HPLN TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)==0) pos_row = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLN_TAN)==0) pos_row = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLN_TAN)==0) pos_row = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLN_TAN)==0) pos_row = 3;

			// HPLT TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLT_TAN)==0) pos_col = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)==0) pos_col = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLT_TAN)==0) pos_col = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLT_TAN)==0) pos_col = 3;			

			// Stokes paramter position , 
			if(strcmp(image->ctype_1,CTYPE_STOKES)==0) pos_stokes_parameters = 0;
			if(strcmp(image->ctype_2,CTYPE_STOKES)==0) pos_stokes_parameters = 1;
			if(strcmp(image->ctype_3,CTYPE_STOKES)==0) pos_stokes_parameters = 2;
			if(strcmp(image->ctype_4,CTYPE_STOKES)==0) pos_stokes_parameters = 3;

			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){


				image->rows=naxes[pos_row];
				image->cols=naxes[pos_col];
				image->nLambdas=naxes[pos_lambda];
				image->numStokes=naxes[pos_stokes_parameters];
				image->numPixels = naxes[pos_col] * naxes[pos_row]; // we will read the image by columns 
				image->pos_lambda = pos_lambda;
				image->pos_col = pos_col;
				image->pos_row = pos_row;
				image->pos_stokes_parameters = pos_stokes_parameters;
				numPixelsFitsFile = naxes[pos_row]*naxes[pos_col]*naxes[pos_lambda]*naxes[pos_stokes_parameters];
				//printf("\n NÚMERO DE PIXELES EN LA IMAGEN %d", numPixelsFitsFile);
				//printf("\n**********************");
				// allocate memory to read all pixels in the same array 
				float * imageTemp = calloc(numPixelsFitsFile, sizeof(float));
				if (!imageTemp)  {
					printf("ERROR ALLOCATION MEMORY FOR TEMP IMAGE");
					return NULL;
          	}
				
				
				long fpixel [4] = {1,1,1,1}; 
				//double time2ReadPixels;
				//clock_t t;
				//t = clock();
				fits_read_pix(fptr, TFLOAT, fpixel, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				//t = clock() - t;
				//time2ReadPixels = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				//printf("\n TIME TO READ PIXELS:  %f seconds to execute \n", time2ReadPixels);				
				if(status){
					fits_report_error(stderr, status);
               return NULL;	
				}

				// allocate memory for reorder the image
				image->pixels = calloc(image->numPixels, sizeof(vpixels));
				if(forParallel){
					//image->vLambdaImagen = calloc(image->numPixels*image->nLambdas, sizeof(PRECISION));
					image->vLambdaImagen = NULL;
					image->spectroImagen = calloc(image->numPixels*image->nLambdas*image->numStokes, sizeof(float));
				}
				else{
					image->vLambdaImagen = NULL;
					image->spectroImagen = NULL;
				}
				//printf("\n Número de pixeles: %d", image->numPixels);
				//printf("\n ***********************************************");
				for( i=0;i<image->numPixels;i++){
					image->pixels[i].spectro = calloc ((image->numStokes*image->nLambdas),sizeof(float));
					//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
					image->pixels[i].vLambda = NULL;
					//image->pixels[i].vLambda = calloc (32, sizeof(PRECISION));
					image->pixels[i].nLambda = image->nLambdas;
				}
				int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0, currentPixel;
				//PRECISION pixel;
				if(naxis==4){ // image with 4 dimension 
					
					for( i=0; i<naxes[3];i++){
						for( j=0; j<naxes[2];j++){
							for( k=0;k<naxes[1];k++){
								for( h=0;h<naxes[0];h++){
									PRECISION pixel = 0.0;
									//fits_read_pix(fptr, datatype, fpixel, 1, &nulval, &pixel, &anynul, &status);
									// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
									switch (pos_lambda)
									{
										case 0:
											currentLambda = h;
											//currentLambda = fpixel[0]-1;
											break;
										case 1:
											currentLambda = k;
											//currentLambda = fpixel[1]-1;
											break;
										case 2:
											currentLambda = j;
											//currentLambda = fpixel[2]-1;
											break;
										case 3:
											currentLambda = i;
											//currentLambda = fpixel[3]-1;
											break;																						
									}
									switch (pos_stokes_parameters)
									{
										case 0:
											currentStokeParameter = h;
											//currentStokeParameter = fpixel[0]-1;
											break;
										case 1:
											currentStokeParameter = k;
											//currentStokeParameter = fpixel[1]-1;
											break;
										case 2:
											currentStokeParameter = j;
											//currentStokeParameter = fpixel[2]-1;
											break;
										case 3:
											currentStokeParameter = i;
											//currentStokeParameter = fpixel[3]-1;
											break;																						
									}
									switch (pos_row)
									{
										case 0:
											currentRow = h;
											//currentRow = fpixel[0]-1;
											break;
										case 1:
											currentRow = k;
											//currentRow = fpixel[1]-1;
											break;
										case 2:
											currentRow = j;
											//currentRow = fpixel[2]-1;
											break;
										case 3:
											currentRow = i;
											//currentRow = fpixel[3]-1;
											break;																						
									}
									switch (pos_col)
									{
										case 0:
											currentCol = h;
											//currentCol = fpixel[0]-1;
											break;
										case 1:
											currentCol = k;
											//currentCol = fpixel[1]-1;;
											break;
										case 2:
											currentCol = j;
											//currentCol = fpixel[2]-1;
											break;
										case 3:
											currentCol = i;
											//currentCol = fpixel[3]-1;
											break;																						
									}			
									pixel = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];
									currentPixel = (currentCol*naxes[pos_row]) + currentRow;
									//currentPixel = (currentRow*naxes[pos_col]) + currentCol;
									//printf("\n CURRENTLAMBDA %d CURRENTSTOKEPARAMETER %d CURRENTROW %d CURRENTCOL %d NUMERO DE SPECTRO %d NÚMERO DE ITER %d -- NUMERO DE PIXEL -- %d  VALOR PIXEL: %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n %d %d %d %d %d %d %d %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n*");
									image->pixels[currentPixel].spectro[currentLambda + (image->nLambdas * currentStokeParameter)] = pixel;  // I =0, Q = 1, U = 2, V = 3
								}
							}
						}
					}
				}
				if(forParallel){
					int contSpectro = 0;
					for( i=0;i<image->numPixels;i++){
						for( j=0;j<(image->nLambdas*image->numStokes);j++){
							image->spectroImagen[contSpectro++] = image->pixels[i].spectro[j];
						}
					}
				}

				free(imageTemp);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return NULL;
				}				
			}
		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); 
      return NULL;
   }
	
	return image; 
}*/

FitsImage *  readFitsSpectroImage (const char * fitsFileSpectra, int forParallel){
   fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
   FitsImage * image =  malloc(sizeof(FitsImage));
   int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, naxis, anynul, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	char comment[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
	
	int i, j, k, h;
   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   if (!fits_open_file(&fptr, fitsFileSpectra, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;

			// ORDER MUST BE CTYPE1->'HPLN-TAN',CTYPE2->'HPLT-TAN',CTYPE3->'WAVE-GRI',CTYPE4->'STOKES'
			// int pos_row = 0, pos_col = 1, pos_lambda = 2, pos_stokes_parameters = 3;
			int correctOrder =0;
			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_lambda; 
			int pos_row;
			int pos_col;
			int pos_stokes_parameters;
			// LAMBDA POSITION
			if(strcmp(image->ctype_1,CTYPE_WAVE)==0) pos_lambda = 0;
			if(strcmp(image->ctype_2,CTYPE_WAVE)==0) pos_lambda = 1;
			if(strcmp(image->ctype_3,CTYPE_WAVE)==0) pos_lambda = 2;
			if(strcmp(image->ctype_4,CTYPE_WAVE)==0) pos_lambda = 3;

			// HPLN TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)==0) pos_row = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLN_TAN)==0) pos_row = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLN_TAN)==0) pos_row = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLN_TAN)==0) pos_row = 3;

			// HPLT TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLT_TAN)==0) pos_col = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)==0) pos_col = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLT_TAN)==0) pos_col = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLT_TAN)==0) pos_col = 3;			

			// Stokes paramter position , 
			if(strcmp(image->ctype_1,CTYPE_STOKES)==0) pos_stokes_parameters = 0;
			if(strcmp(image->ctype_2,CTYPE_STOKES)==0) pos_stokes_parameters = 1;
			if(strcmp(image->ctype_3,CTYPE_STOKES)==0) pos_stokes_parameters = 2;
			if(strcmp(image->ctype_4,CTYPE_STOKES)==0) pos_stokes_parameters = 3;

			
			
			if(pos_row==0 && pos_col ==1 && pos_lambda == 2 && pos_stokes_parameters == 3)
			//if(pos_row==2 && pos_col ==3 && pos_lambda == 0 && pos_stokes_parameters == 1)
				correctOrder = 1;
			// GET THE CURRENT POSITION OF EVERY PARAMETER
			/*int pos_row = 2;
			int pos_col = 3;
			int pos_lambda = 0; 
			int pos_stokes_parameters = 1;*/

			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){

				if(bitpix != FLOAT_IMG){
					printf("\n ERROR: the datatype of FITS spectro image must be FLOAT\n");
					printf("\n EXITING THE PROGRAM");
					fits_close_file(fptr, &status);
					exit(EXIT_FAILURE);
				}



				image->rows=naxes[pos_row];
				image->cols=naxes[pos_col];
				image->nLambdas=naxes[pos_lambda];
				image->numStokes=naxes[pos_stokes_parameters];
				image->numPixels = naxes[pos_col] * naxes[pos_row]; // we will read the image by columns 
				image->pos_lambda = pos_lambda;
				image->pos_col = pos_col;
				image->pos_row = pos_row;
				image->pos_stokes_parameters = pos_stokes_parameters;
				numPixelsFitsFile = naxes[pos_row]*naxes[pos_col]*naxes[pos_lambda]*naxes[pos_stokes_parameters];

				// allocate memory to read all pixels in the same array 
				float * imageTemp = calloc(numPixelsFitsFile, sizeof(float));
				if (!imageTemp)  {
					printf("ERROR ALLOCATION MEMORY FOR TEMP IMAGE");
					return NULL;
          	}
				
				
				long fpixel [4] = {1,1,1,1}; 
				fits_read_pix(fptr, TFLOAT, fpixel, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				if(status){
					fits_report_error(stderr, status);
               return NULL;	
				}

				// allocate memory for reorder the image
				
				if(forParallel){
					//image->vLambdaImagen = calloc(image->numPixels*image->nLambdas, sizeof(PRECISION));
					image->vLambdaImagen = NULL;
					image->pixels = NULL;
					image->spectroImagen = calloc(image->numPixels*image->nLambdas*image->numStokes, sizeof(float));
				}
				else{
					image->pixels = calloc(image->numPixels, sizeof(vpixels));
					for( i=0;i<image->numPixels;i++){
						image->pixels[i].spectro = calloc ((image->numStokes*image->nLambdas),sizeof(float));
						//image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
						image->pixels[i].vLambda = NULL;
						image->pixels[i].nLambda = image->nLambdas;
					}					
					image->vLambdaImagen = NULL;
					image->spectroImagen = NULL;
				}


				
				//PRECISION pixel;
				if(naxis==4){ // image with 4 dimension 
					if(correctOrder){
						// i = stokes, j = lambda, k = cols, h = rows 		
						/*for( i=0; i<naxes[3];i++){
							for( j=0; j<naxes[2];j++){
								for( k=0;k<naxes[1];k++){
									for( h=0;h<naxes[0];h++){
										image->pixels[(k*image->rows) + h].spectro[j + (image->nLambdas * i)] = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];  // I =0, Q = 1, U = 2, V = 3
									}
								}
							}
						}*/
						// i = cols, j = rows, k = stokes, h = lambda
						for( i=0; i<naxes[3];i++){
							for( j=0; j<naxes[2];j++){
								for( k=0;k<naxes[1];k++){
									for( h=0;h<naxes[0];h++){
										if(forParallel){
											image->spectroImagen[ (((i*naxes[2]) + j)*(image->nLambdas*image->numStokes)) + (image->nLambdas * k) + h] = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];  // I =0, Q = 1, U = 2, V = 3
										}
										else{
											image->pixels[(i*naxes[2]) + j].spectro[h + (image->nLambdas * k)] = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];  // I =0, Q = 1, U = 2, V = 3
										}
									}
								}
							}
						}					
					}
					else{
						int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0, currentPixel;
						for( i=0; i<naxes[3];i++){
							for( j=0; j<naxes[2];j++){
								for( k=0;k<naxes[1];k++){
									for( h=0;h<naxes[0];h++){
										PRECISION pixel = 0.0;
										//fits_read_pix(fptr, datatype, fpixel, 1, &nulval, &pixel, &anynul, &status);
										// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
										switch (pos_lambda)
										{
											case 0:
												currentLambda = h;
												//currentLambda = fpixel[0]-1;
												break;
											case 1:
												currentLambda = k;
												//currentLambda = fpixel[1]-1;
												break;
											case 2:
												currentLambda = j;
												//currentLambda = fpixel[2]-1;
												break;
											case 3:
												currentLambda = i;
												//currentLambda = fpixel[3]-1;
												break;																						
										}
										switch (pos_stokes_parameters)
										{
											case 0:
												currentStokeParameter = h;
												//currentStokeParameter = fpixel[0]-1;
												break;
											case 1:
												currentStokeParameter = k;
												//currentStokeParameter = fpixel[1]-1;
												break;
											case 2:
												currentStokeParameter = j;
												//currentStokeParameter = fpixel[2]-1;
												break;
											case 3:
												currentStokeParameter = i;
												//currentStokeParameter = fpixel[3]-1;
												break;																						
										}
										switch (pos_row)
										{
											case 0:
												currentRow = h;
												//currentRow = fpixel[0]-1;
												break;
											case 1:
												currentRow = k;
												//currentRow = fpixel[1]-1;
												break;
											case 2:
												currentRow = j;
												//currentRow = fpixel[2]-1;
												break;
											case 3:
												currentRow = i;
												//currentRow = fpixel[3]-1;
												break;																						
										}
										switch (pos_col)
										{
											case 0:
												currentCol = h;
												//currentCol = fpixel[0]-1;
												break;
											case 1:
												currentCol = k;
												//currentCol = fpixel[1]-1;;
												break;
											case 2:
												currentCol = j;
												//currentCol = fpixel[2]-1;
												break;
											case 3:
												currentCol = i;
												//currentCol = fpixel[3]-1;
												break;																						
										}			
										pixel = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];
										currentPixel = (currentCol*naxes[pos_row]) + currentRow;
										if(forParallel){
											image->spectroImagen[(currentPixel *(image->nLambdas*image->numStokes)) + (image->nLambdas * currentStokeParameter) + currentLambda] = pixel;
										}
										else
										{
											image->pixels[currentPixel].spectro[currentLambda + (image->nLambdas * currentStokeParameter)] = pixel;  // I =0, Q = 1, U = 2, V = 3
										}
										
										
									}
								}
							}
						}

					}
				}
				free(imageTemp);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return NULL;
				}				
			}
		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); /* print any error message */
      return NULL;
   }
	
	return image; 
}



/**
 * @param fitsFileSpectra: name of file with FITS image. 
 * @param configControlFile: pointer to control structure with information relative to whole configuration of program. 
 * 
 * */

FitsImage * readFitsSpectroNPixels (const char * fitsFileSpectra, int rowInit, int colInit,int rowEnd, int colEnd, int numPixelsRead,  int forParallel){



	fitsfile *fptr;   
	FitsImage * image =  malloc(sizeof(FitsImage));
   int status = 0;   
	PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, naxis, anynul, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; 
	char comment[FLEN_CARD];
	
	int i, j, k, h;
   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   if (!fits_open_file(&fptr, fitsFileSpectra, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;


			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_lambda; 
			int pos_row;
			int pos_col;
			int pos_stokes_parameters;
			// LAMBDA POSITION
			if(strcmp(image->ctype_1,CTYPE_WAVE)==0) pos_lambda = 0;
			if(strcmp(image->ctype_2,CTYPE_WAVE)==0) pos_lambda = 1;
			if(strcmp(image->ctype_3,CTYPE_WAVE)==0) pos_lambda = 2;
			if(strcmp(image->ctype_4,CTYPE_WAVE)==0) pos_lambda = 3;

			// HPLN TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)==0) pos_row = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLN_TAN)==0) pos_row = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLN_TAN)==0) pos_row = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLN_TAN)==0) pos_row = 3;

			// HPLT TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLT_TAN)==0) pos_col = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)==0) pos_col = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLT_TAN)==0) pos_col = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLT_TAN)==0) pos_col = 3;			

			// Stokes paramter position , 
			if(strcmp(image->ctype_1,CTYPE_STOKES)==0) pos_stokes_parameters = 0;
			if(strcmp(image->ctype_2,CTYPE_STOKES)==0) pos_stokes_parameters = 1;
			if(strcmp(image->ctype_3,CTYPE_STOKES)==0) pos_stokes_parameters = 2;
			if(strcmp(image->ctype_4,CTYPE_STOKES)==0) pos_stokes_parameters = 3;

			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){


				//image->rows=naxes[pos_row];
				//image->rows=(subx2-subx1)+1;
				//image->cols=naxes[pos_col];
				//image->cols= (suby2-suby1)+1;
				image->nLambdas=naxes[pos_lambda];
				image->numStokes=naxes[pos_stokes_parameters];
				//image->numPixels = naxes[pos_col] * naxes[pos_row]; // we will read the image by columns 
				image->numPixels = numPixelsRead ; // we will read the image by columns 
				image->pos_lambda = pos_lambda;
				image->pos_col = pos_col;
				image->pos_row = pos_row;
				image->pos_stokes_parameters = pos_stokes_parameters;
				numPixelsFitsFile = image->numPixels*image->nLambdas*image->numStokes;
				//printf("\n NÚMERO DE PIXELES EN LA IMAGEN %d", numPixelsFitsFile);
				//printf("\n**********************");
				// allocate memory to read all pixels in the same array 
				float * imageTemp = calloc(numPixelsFitsFile, sizeof(float));
				if (!imageTemp)  {
					printf("ERROR ALLOCATION MEMORY FOR TEMP IMAGE");
					return NULL;
          	}
				
				
				long fpixelBegin [4] = {1,1,1,1}; 
				//long fpixelEnd [4] = {1,1,1,1}; 
				//long inc [4] = {1,1,1,1};
				fpixelBegin[pos_row] = rowInit;
				//fpixelEnd[pos_row] = subx2;

				fpixelBegin[pos_col] = colInit;
				//fpixelEnd[pos_col] = suby2;

				//fpixelEnd[pos_lambda] = naxes[pos_lambda];
				//fpixelEnd[pos_stokes_parameters] = naxes[pos_stokes_parameters];


				double time2ReadPixels;
				clock_t t;
				t = clock();
				
				fits_read_pix(fptr, TFLOAT, fpixelBegin, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				t = clock() - t;
				time2ReadPixels = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n ID PROC: %d TIME TO READ PIXELS:  %f seconds to execute \n",forParallel, time2ReadPixels);				
				if(status){
					fits_report_error(stderr, status);
            		return NULL;	
				}
				
				// allocate memory for reorder the image
				
				if(forParallel){
					//image->vLambdaImagen = calloc(image->numPixels*image->nLambdas, sizeof(PRECISION));
					image->vLambdaImagen = NULL;
					image->spectroImagen = calloc(image->numPixels*image->nLambdas*image->numStokes, sizeof(float));
				}
				else{
					image->vLambdaImagen = NULL;
					image->spectroImagen = NULL;
					image->pixels = calloc(image->numPixels, sizeof(vpixels));
					for( i=0;i<image->numPixels;i++){
						image->pixels[i].spectro = calloc ((image->numStokes*image->nLambdas),sizeof(float));
						image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
						//image->pixels[i].vLambda = calloc (32, sizeof(PRECISION));
						image->pixels[i].nLambda = image->nLambdas;
					}					
				}
				printf("\n IDPROD %d Número de pixeles: %d",forParallel, image->numPixels);
				printf("\n ***********************************************");

				int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0, currentPixel;
				//PRECISION pixel;
				if(naxis==4){ // image with 4 dimension 
					
					for( i=0; i<naxes[3];i++){
						for( j=0; j<naxes[2];j++){
							for( k=0;k<naxes[1];k++){
								for( h=0;h<naxes[0];h++){
//					for( fpixel[3] = 1; fpixel[3]<=naxes[3];fpixel[3]++){
//						for( fpixel[2] = 1; fpixel[2]<=naxes[2];fpixel[2]++){
//							for( fpixel[1] = 1; fpixel[1]<=naxes[1];fpixel[1]++){
//								for( fpixel[0] = 1; fpixel[0]<=naxes[0]; fpixel[0]++){	
									PRECISION pixel = 0.0;
									//fits_read_pix(fptr, datatype, fpixel, 1, &nulval, &pixel, &anynul, &status);
									// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
									switch (pos_lambda)
									{
										case 0:
											currentLambda = h;
											//currentLambda = fpixel[0]-1;
											break;
										case 1:
											currentLambda = k;
											//currentLambda = fpixel[1]-1;
											break;
										case 2:
											currentLambda = j;
											//currentLambda = fpixel[2]-1;
											break;
										case 3:
											currentLambda = i;
											//currentLambda = fpixel[3]-1;
											break;																						
									}
									switch (pos_stokes_parameters)
									{
										case 0:
											currentStokeParameter = h;
											//currentStokeParameter = fpixel[0]-1;
											break;
										case 1:
											currentStokeParameter = k;
											//currentStokeParameter = fpixel[1]-1;
											break;
										case 2:
											currentStokeParameter = j;
											//currentStokeParameter = fpixel[2]-1;
											break;
										case 3:
											currentStokeParameter = i;
											//currentStokeParameter = fpixel[3]-1;
											break;																						
									}
									switch (pos_row)
									{
										case 0:
											currentRow = h;
											//currentRow = fpixel[0]-1;
											break;
										case 1:
											currentRow = k;
											//currentRow = fpixel[1]-1;
											break;
										case 2:
											currentRow = j;
											//currentRow = fpixel[2]-1;
											break;
										case 3:
											currentRow = i;
											//currentRow = fpixel[3]-1;
											break;																						
									}
									switch (pos_col)
									{
										case 0:
											currentCol = h;
											//currentCol = fpixel[0]-1;
											break;
										case 1:
											currentCol = k;
											//currentCol = fpixel[1]-1;;
											break;
										case 2:
											currentCol = j;
											//currentCol = fpixel[2]-1;
											break;
										case 3:
											currentCol = i;
											//currentCol = fpixel[3]-1;
											break;																						
									}			
									pixel = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];
									currentPixel = (currentCol*naxes[pos_row]) + currentRow;
									if(forParallel){
										image->spectroImagen[currentPixel*(image->nLambdas*image->numStokes)+(image->nLambdas * currentStokeParameter)+currentLambda] = pixel;
									}
									else
									{
										image->pixels[currentPixel].spectro[currentLambda + (image->nLambdas * currentStokeParameter)] = pixel;  // I =0, Q = 1, U = 2, V = 3
									}
									
									//currentPixel = (currentRow*naxes[pos_col]) + currentCol;
									//printf("\n CURRENTLAMBDA %d CURRENTSTOKEPARAMETER %d CURRENTROW %d CURRENTCOL %d NUMERO DE SPECTRO %d NÚMERO DE ITER %d -- NUMERO DE PIXEL -- %d  VALOR PIXEL: %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n %d %d %d %d %d %d %d %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n*");
									
								}
							}
						}
					}
				}

				free(imageTemp);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return NULL;
				}				
			}
		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); 
      return NULL;
   }
	
	return image; 
}

/*FitsImage * readFitsSpectroImageRectangular (const char * fitsFileSpectra, ConfigControl * configCrontrolFile, int forParallel){



	fitsfile *fptr;   
   FitsImage * image =  malloc(sizeof(FitsImage));
   int status = 0;   
   PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, naxis, anynul, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; 
	char comment[FLEN_CARD];   
	
	int i, j, k, h;
   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   if (!fits_open_file(&fptr, fitsFileSpectra, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;

			// ORDER MUST BE CTYPE1->'HPLN-TAN',CTYPE2->'HPLT-TAN',CTYPE3->'WAVE-GRI',CTYPE4->'STOKES'

			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)!=0){
				printf("\n ERROR reading spectro image, the first dimension must be 'HPLN-TAN' and is %s",image->ctype_1);
			}
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)!=0){
				printf("\n ERROR reading spectro image, the first dimension must be 'HPLN-TAN' and is %s",image->ctype_1);
			}
			if(strcmp(image->ctype_3,CTYPE_WAVE)!=0){
				printf("\n ERROR reading spectro image, the first dimension must be 'HPLN-TAN' and is %s",image->ctype_1);
			}
			if(strcmp(image->ctype_4,CTYPE_STOKES)!=0){
				printf("\n ERROR reading spectro image, the first dimension must be 'HPLN-TAN' and is %s",image->ctype_1);
			}

			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_row = 0;
			int pos_col = 1;
			int pos_lambda = 2; 
			int pos_stokes_parameters = 3;


			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){


				//image->rows=naxes[pos_row];
				if(configCrontrolFile->subx2==0 && configCrontrolFile->subx1==0){
					image->rows = naxes[0];
					configCrontrolFile->subx1 = 1;
					configCrontrolFile->subx2 = naxes[0];
				}
				else			
					image->rows=(configCrontrolFile->subx2-configCrontrolFile->subx1)+1;
				
				if(configCrontrolFile->suby2==0 && configCrontrolFile->suby1==0){
					image->cols = naxes[1];
					configCrontrolFile->suby1 = 1;
					configCrontrolFile->suby2 = naxes[1];
				}
				else
					image->cols= (configCrontrolFile->suby2-configCrontrolFile->suby1)+1;

				image->nLambdas=naxes[2];
				image->numStokes=naxes[3];

				image->numPixels = image->cols * image->rows ; // we will read the image by columns 
				image->pos_row = 0;
				image->pos_col = 1;
				image->pos_lambda = 2;
				image->pos_stokes_parameters = 3;
				numPixelsFitsFile = image->rows*image->rows*image->nLambdas*image->numStokes;
				
				
				// allocate memory to read all pixels in the same array 
				float * imageTemp = calloc(numPixelsFitsFile, sizeof(float));
				if (!imageTemp)  {
					printf("ERROR ALLOCATION MEMORY FOR TEMP IMAGE");
					return NULL;
          	}
				
				
				long fpixelBegin [4] = {1,1,1,1}; 
				long fpixelEnd [4] = {1,1,1,1}; 
				long inc [4] = {1,1,1,1};
				fpixelBegin[pos_row] = configCrontrolFile->subx1;
				fpixelEnd[pos_row] = configCrontrolFile->subx2;

				fpixelBegin[pos_col] = configCrontrolFile->suby1;
				fpixelEnd[pos_col] = configCrontrolFile->suby2;

				fpixelEnd[pos_lambda] = naxes[pos_lambda];
				fpixelEnd[pos_stokes_parameters] = naxes[pos_stokes_parameters];


				//double time2ReadPixels;
				//clock_t t;
				//t = clock();
				fits_read_subset(fptr, TFLOAT, fpixelBegin, fpixelEnd, inc, &nulval, imageTemp, &anynul, &status);
				//fits_read_pix(fptr, TDOUBLE, fpixel, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				//t = clock() - t;
				//time2ReadPixels = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				//printf("\n TIME TO READ PIXELS:  %f seconds to execute \n", time2ReadPixels);				
				if(status){
					fits_report_error(stderr, status);
               return NULL;	
				}

				// allocate memory for reorder the image
				image->pixels = calloc(image->numPixels, sizeof(vpixels));
				if(forParallel){
					image->vLambdaImagen = calloc(image->numPixels*image->nLambdas, sizeof(float));
					image->spectroImagen = calloc(image->numPixels*image->nLambdas*image->numStokes, sizeof(float));
				}
				else{
					image->vLambdaImagen = NULL;
					image->spectroImagen = NULL;
				}
				//printf("\n Número de pixeles: %d", image->numPixels);
				//printf("\n ***********************************************");
				for( i=0;i<image->numPixels;i++){
					image->pixels[i].spectro = calloc ((image->numStokes*image->nLambdas),sizeof(float));
					image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(float));
					image->pixels[i].nLambda = image->nLambdas;
				}
				
				//PRECISION pixel;
				if(naxis==4){ // image with 4 dimension 
					
					int sizeDim0 = (fpixelEnd[0]-(fpixelBegin[0]-1));
					int sizeDim1 = (fpixelEnd[1]-(fpixelBegin[1]-1));
					int sizeDim2 = (fpixelEnd[2]-(fpixelBegin[2]-1));
					int sizeDim3 = (fpixelEnd[3]-(fpixelBegin[3]-1));
					

					// i = stokes, j = lambda, k = cols, h = rows 
					for( i=0; i<sizeDim3;i++){
						for( j=0; j<sizeDim2;j++){
							for( k=0;k<sizeDim1;k++){
								for( h=0;h<sizeDim0;h++){
									//currentPixel = (k*image->rows) + h;
									image->pixels[(k*image->rows) + h].spectro[j + (image->nLambdas * i)] = imageTemp [(i*(fpixelEnd[2]-(fpixelBegin[2]-1))*(fpixelEnd[1]-(fpixelBegin[1]-1))*(fpixelEnd[0]-(fpixelBegin[0]-1))) + (j*(fpixelEnd[1]-(fpixelBegin[1]-1))*(fpixelEnd[0]-(fpixelBegin[0]-1))) + (k*(fpixelEnd[0]-(fpixelBegin[0]-1))) + h];
								}
							}
						}
					}
				}
				if(forParallel){
					int contSpectro = 0;
					for( i=0;i<image->numPixels;i++){
						for( j=0;j<(image->nLambdas*image->numStokes);j++){
							image->spectroImagen[contSpectro++] = image->pixels[i].spectro[j];
						}
					}
				}

				free(imageTemp);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return NULL;
				}				
			}
		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); 
      return NULL;
   }
	
	return image; 
}
*/



/**
 * This function read the lambda values for the image from the file "fitsFileLambda" and store it into a struct of FitsImage. The file of spectro must
 * be read it before call this method. 
 * fitsFileLambda --> name of the fits file to read with lambda values 
 * fitsImage --> struct of image 
 * Return 1 If the image has been read corectly if not return 0 
 */

int  readFitsLambdaFile (const char * fitsFileLambda, FitsImage * fitsImage){
	int i, j, k;
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION  nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
	int bitpix, naxis, anynul;
	long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	
	/*printf("\n READING IMAGE WITH LAMBDA ");
	printf("\n**********");*/
	if (!fits_open_file(&fptr, fitsFileLambda, READONLY, &status)){
		if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
			/*int datatype = 0;
			switch(bitpix) {
				case BYTE_IMG:
					datatype = TBYTE;
					break;
				case SHORT_IMG:
					datatype = TSHORT;
					break;
				case LONG_IMG:
					datatype = TINT;
					break;
				case FLOAT_IMG:
					datatype = TFLOAT;
					break;
				case DOUBLE_IMG:
					datatype = TDOUBLE;
					break;
			}*/

			if(naxis!=1  || naxis!=3){
				if(naxis == 1){ // array of lambads 
					if(naxes[0]!=fitsImage->nLambdas){ // image of lambas has different size of spectra image 
						printf("\n IMAGE OF LAMBAS HAS DIFFERENT SIZE OF SPECTRA IMAGE. NUMBER OF LAMBDAS IN SPECTRA IMAGE %d NUMBER OF LAMBDA IMAGE %ld ", fitsImage->nLambdas,naxes[0]);
						freeFitsImage(fitsImage);
						return 0;
					}
					PRECISION  * vAuxLambdas = calloc(naxes[0], sizeof(PRECISION));
					long fpixel [1] = {1};
					fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], &nulval, vAuxLambdas, &anynul, &status) ;
					if(status){
						fits_report_error(stderr, status);
						return 0;	
					}
					int contLambda = 0;

					for( i=0;i<fitsImage->numPixels;i++){
						for( j=0;j<naxes[0];j++){
							fitsImage->pixels[i].vLambda[j]=vAuxLambdas[j];
							if(fitsImage->vLambdaImagen!=NULL)
								fitsImage->vLambdaImagen[contLambda++] = fitsImage->pixels[i].vLambda[j];
						}
					}								
				}
				else if(naxis == 3){  // matrix of lambdas  
					if( naxes[0]!= fitsImage->rows || naxes[1]!= fitsImage->cols || naxes[2]!=fitsImage->nLambdas){ // image of lambas has different size of spectra image 
						printf("\n IMAGE OF LAMBAS HAS DIFFERENT SIZE OF SPECTRA IMAGE. SIZE SPECTRA %d X %d X %d. SIZE LAMBDA IMAGE %ld X %ld X %ld", fitsImage->rows, fitsImage->cols, fitsImage->nLambdas,naxes[0], naxes[1], naxes[2]);
						freeFitsImage(fitsImage);
						return 0;
					}
					// READ ALL FILE IN ONLY ONE ARRAY 
					// WE ASSUME THAT DATA COMES IN THE FORMAT ROW x COL x LAMBDA
					int numLambdas2Read = naxes[0]*naxes[1]*naxes[2];
					PRECISION  * vAuxLambdas = calloc(numLambdas2Read, sizeof(PRECISION));
					
					//fits_read_img(fptr, datatype, first, numLambdas2Read, &nulval, vAuxLambdas, &anynul, &status);
					long fpixel [3] = {1,1,1};
					fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1]*naxes[2], &nulval, vAuxLambdas, &anynul, &status);
					if(status){
						fits_report_error(stderr, status);
						return 0;	
					}
					int offset_1 = naxes[1] * naxes[0];
					int offset_2 = naxes[0];
					int contLambda = 0;
					for( i=0;i<naxes[2];i++){ // LAMBDA
						for( j=0;j<naxes[1];j++){ // COLS 
							for( k=0;j<naxes[0];k++){ // ROWS
								fitsImage->pixels[ ((j* offset_2) + k) ].vLambda[i] = vAuxLambdas [ (i*offset_1) + (j*offset_2)+ k];
								fitsImage->vLambdaImagen[contLambda++] = fitsImage->pixels[ ((j* offset_2) + k) ].vLambda[i];
							}
						}
					}
				}
			}
			else{
				printf("\n NAXIS FROM LAMBA FILE IS NOT VALID %d ** \n", naxis);
				freeFitsImage(fitsImage);
				return 0;
			}
			// CLOSE FILE FITS LAMBDAS
			fits_close_file(fptr, &status);
			if (status){
				fits_report_error(stderr, status);
				return 0;
			}
		}
		else {
			printf("\n WE CAN NOT OPEN FILE OF LAMBAS ** \n");
			if (status) fits_report_error(stderr, status); /* print any error message */
			freeFitsImage(fitsImage);
			return 0;
		}
	}
	else {
		printf("\n WE CAN NOT READ PARAMETERS FROM THE FILE  %s \n",fitsFileLambda);
		if (status) fits_report_error(stderr, status); /* print any error message */
		freeFitsImage(fitsImage);
		return 0;
	}
	/*printf("\n LAMBDA IMAGE READ");
	printf("\n**********");*/

	return 1;

}

/**
 * This function read the lambda values for the image from the file "fitsFileLambda" and store it into a struct of FitsImage. The file of spectro must
 * be read it before call this method. 
 * fitsFileLambda --> name of the fits file to read with lambda values 
 * fitsImage --> struct of image 
 * @return Vector with 
 */

PRECISION * readFitsLambdaToArray (const char * fitsFileLambda,  int * indexLine, int * nLambda){
	int i;
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION  nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
	PRECISION * vLambda = NULL;
	int bitpix, naxis, anynul;
	long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	
	/*printf("\n READING IMAGE WITH LAMBDA ");
	printf("\n**********");*/
	if (!fits_open_file(&fptr, fitsFileLambda, READONLY, &status)){
		if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
			if(naxis!=2  || naxis!=4){
				if(naxis == 2){ // array of lambads 

					*nLambda = naxes[0];
					long fpixel [2] = {1,1};
					i=0;
					vLambda = calloc(*nLambda,sizeof(PRECISION));
					for(fpixel[1]=1;fpixel[1]<=naxes[1];fpixel[1]++){
						for(fpixel[0]=1;fpixel[0]<=naxes[0];fpixel[0]++){
							PRECISION lambdaAux;
							fits_read_pix(fptr, TDOUBLE, fpixel, 1, &nulval, &lambdaAux, &anynul, &status) ;		
							if(fpixel[1]==1)
								*indexLine = (int) lambdaAux;
							else
								vLambda[i++] = lambdaAux;
						}
					}
					printf("\nLAMBDAS LEIDOS:\n");
					for(i=0;i<*nLambda;i++){
						printf(" %lf ",vLambda[i]);
					}
					printf("\n");
					//fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], &nulval, data, &anynul, &status) ;
					if(status){
						fits_report_error(stderr, status);
						free(vLambda);
						return 0;	
					}
				}
				else if(naxis == 4){  // matrix of lambdas  
					// READ ALL FILE IN ONLY ONE ARRAY 
					// WE ASSUME THAT DATA COMES IN THE FORMAT ROW x COL x LAMBDA					
					int numLambdas2Read = naxes[0]*naxes[1]*naxes[2];
					//fits_read_img(fptr, datatype, first, numLambdas2Read, &nulval, vAuxLambdas, &anynul, &status);
					long fpixel [3] = {1,1,1};
					fits_read_pix(fptr, TDOUBLE, fpixel, numLambdas2Read, &nulval, vLambda, &anynul, &status);
					if(status){
						fits_report_error(stderr, status);
						free(vLambda);
						return 0;	
					}
				}
			}
			else{
				printf("\n NAXIS FROM LAMBA FILE IS NOT VALID %d ** \n", naxis);
				return 0;
			}
			// CLOSE FILE FITS LAMBDAS
			fits_close_file(fptr, &status);
			if (status){
				fits_report_error(stderr, status);
				return 0;
			}
		}
		else {
			printf("\n WE CAN NOT OPEN FILE OF LAMBAS ** \n");
			if (status) fits_report_error(stderr, status); /* print any error message */
			return 0;
		}
	}
	else {
		printf("\n WE CAN NOT READ PARAMETERS FROM THE FILE  %s \n",fitsFileLambda);
		if (status) fits_report_error(stderr, status); /* print any error message */
		return 0;
	}
	/*printf("\n LAMBDA IMAGE READ");
	printf("\n**********");*/

	return vLambda;

}



PRECISION * readFitsStrayLightFile (const char * fitsFileStrayLight, int * dimStrayLight, int numLambda){
	
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
	int bitpix, naxis, anynul;
	long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	
	PRECISION * vStrayLight = NULL;
	/*printf("\n READING IMAGE WITH LAMBDA ");
	printf("\n**********");*/
	if (!fits_open_file(&fptr, fitsFileStrayLight, READONLY, &status)){
		if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
			/*int datatype = 0;
			switch(bitpix) {
				case BYTE_IMG:
					datatype = TBYTE;
					break;
				case SHORT_IMG:
					datatype = TSHORT;
					break;
				case LONG_IMG:
					datatype = TINT;
					break;
				case FLOAT_IMG:
					datatype = TFLOAT;
					break;
				case DOUBLE_IMG:
					datatype = TDOUBLE;
					break;
			}*/

			if(naxis==2){
				if( naxes[0]!= numLambda || naxes[1]!= NPARMS){ // stray light has different size of spectra image 
					printf("\n STRAY LIGHT FILE HAS DIFFERENT SIZE OF SPECTRA IMAGE. SIZE SPECTRA %d X %d X . STRAY LIGHT SIZE %ld X %ld ", numLambda, NPARMS , naxes[0], naxes[1]);
					return NULL;
				}
				// READ ALL FILE IN ONLY ONE ARRAY 
				// WE ASSUME THAT DATA COMES IN THE FORMAT ROW x COL x LAMBDA
				*dimStrayLight = naxes[0]*naxes[1];
				vStrayLight = calloc(*dimStrayLight, sizeof(PRECISION));
				long fpixel [3] = {1,1,1};
				fits_read_pix(fptr, TDOUBLE, fpixel, *dimStrayLight, &nulval, vStrayLight, &anynul, &status);
				if(status){
					fits_report_error(stderr, status);
					return NULL;	
				}
				printf("\n STRAY LIGHT LEIDO: \n");
			}
			else{
				printf("\n NAXIS FROM STRAY LIGHT FILE IS NOT VALID %d ** \n", naxis);
				return 0;
			}
			// CLOSE FILE FITS LAMBDAS
			fits_close_file(fptr, &status);
			if (status){
				fits_report_error(stderr, status);
				return  NULL;
			}
		}
		else {
			printf("\n WE CAN NOT OPEN FILE OF STRAY LIGHT ** \n");
			if (status) fits_report_error(stderr, status); /* print any error message */
			return NULL;
		}
	}
	else {
		printf("\n WE CAN NOT READ PARAMETERS FROM THE FILE  %s \n",fitsFileStrayLight);
		if (status) fits_report_error(stderr, status); /* print any error message */
		return NULL;
	}
	printf("\n STRAY LIGHT FILE READ");
	printf("\n**********");

	return vStrayLight;
}


void freeFitsImage(FitsImage * image){
	int i;
	if(image!=NULL){
		if(image->pixels!=NULL){
			for( i=0;i<image->numPixels;i++){
				if(image->pixels[i].spectro!=NULL)
					free(image->pixels[i].spectro);
				if(image->pixels[i].vLambda!=NULL)
					free(image->pixels[i].vLambda);
			}

			free(image->pixels);
		}
		if(image->spectroImagen!=NULL)
			free(image->spectroImagen);
		if(image->vLambdaImagen!=NULL)
			free(image->vLambdaImagen);
		free(image);
	}
	
}


/**
 * 
 * fixed = array with positions to write in the file, Positions are in the following order: 
 * [Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
 * */
int writeFitsImageModels(const char * fitsFile, int numRows, int numCols, Init_Model * vInitModel, float * vChisqrf, int * vNumIterPixel, int addChiqr){

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status;
	int i, j, h; // indexes for loops
	long  fpixel;
	int indexModel = 0; 


	int bitpix =  FLOAT_IMG; 
	long naxis =   3;  /* 2-dimensional image */    
	long naxes[3] = { numRows, numCols, NUMBER_PARAM_MODELS+1};   /* Image of numRows X numCols x 10 parameters of model and chisqrf */


	if(addChiqr){
		naxes[2]++;
	}
   

   remove(fitsFile);               /* Delete old file if it already exists */
   status = 0;         /* initialize status before calling fitsio routines */
   if (fits_create_file(&fptr, fitsFile, &status)) /* create new FITS file */
   	printerror( status );           /* call printerror if error occurs */
	
	 /* write the required keywords for the primary array image.     */
    /* Since bitpix = FLOAT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = -32 (float) .Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
		printerror( status );
		return 0;
	}

	float * vModel = calloc(naxes[0] * naxes[1] * naxes[2], sizeof(float));

	for( i=0;i<naxes[2];i++){
		for( j=0;j<naxes[0];j++){
			for( h=0; h<naxes[1];h++){
				//[Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
				switch (i)
				{
				case 0:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].eta0;
					break;
				case 1:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].B;
					break;
				case 2:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].vlos;
					break;
				case 3:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].dopp;
					break;
				case 4:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].aa;
					break;
				case 5:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].gm;
					break;					
				case 6:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].az;
					break;					
				case 7:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].S0;
					break;					
				case 8:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].S1;
					break;					
				case 9:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].mac;
					break;					
				case 10:
					vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].alfa;
					break;
				case 11: // NUMBER OF ITERATIONS
					vModel[indexModel++] = vNumIterPixel[( j*naxes[1]) + h];
					break;
				case 12: // CHISQR 
					vModel[indexModel++] = vChisqrf[( j*naxes[1]) + h];
					break;										
				default:
					break;
				}
			}
		}
	}

   fpixel = 1;                               /* first pixel to write      */
   //nelements = naxes[0] * naxes[1] * naxes[2];          /* number of pixels to write */


   if ( fits_write_img(fptr, TFLOAT, fpixel, indexModel, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}

	// CLEAN MEMORY 
	free(vModel);

	    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
   /*exposure = 1500;
	if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
		"Total Exposure Time", &status) ){
		printerror( status );           
		return 0;
	}*/
	
	if ( fits_close_file(fptr, &status) ){        
		printerror( status );
		return 0;
	}
	
	return 1;

}

/**
 * 
 * fixed = array with positions to write in the file, Positions are in the following order: 
 * [Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
 * */
int writeFitsImageModelsSubSet(const char * fitsFileName, int numRows, int numCols,int rowInit, int colInit, int numPixelsWrite, int displs, Init_Model * vInitModel, float * vChisqrf, int * vNumIterPixel, int addChiqr){
	
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status;
	int i, j, h; // indexes for loops
	
	int indexModel = 0; 
	int bitpix =  FLOAT_IMG; 
	long naxis =   3;  /* 2-dimensional image */    
	long naxes[3] = { numRows, numCols, NUMBER_PARAM_MODELS+1};   /* Image of numRows X numCols x 10 parameters of model and chisqrf */
	if(addChiqr){
		naxes[2]++;
	}	
	// check if the file not exist and create the image, close inmediatly to prevent another process to write
	if(access(fitsFileName,F_OK) == -1){
		remove(fitsFileName);               // Delete old file if it already exists 
		status = 0;         // initialize status before calling fitsio routines 
		if (fits_create_file(&fptr, fitsFileName, &status)) // create new FITS file 
			printerror( status );           // call printerror if error occurs 
		

		if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
			printerror( status );
			return 0;
		}
		if ( fits_close_file(fptr, &status) ){        
			printerror( status );
			return 0;
		}
		printf("\nfichero creado de nuevo y cerrado");
	}

	// open the file always like image
	fits_open_file(&fptr, fitsFileName, READWRITE, &status);

	
	float * vModel = calloc(numPixelsWrite * naxes[2], sizeof(float));

	for(i=0;i<naxes[2];i++){
		for(j=0;j<numPixelsWrite;j++){
			//[Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
			switch (i)
			{
			case 0:
				vModel[indexModel++] = vInitModel[j].eta0;
				break;
			case 1:
				vModel[indexModel++] = vInitModel[j].B;
				break;
			case 2:
				vModel[indexModel++] = vInitModel[j].vlos;
				break;
			case 3:
				vModel[indexModel++] = vInitModel[j].dopp;
				break;
			case 4:
				vModel[indexModel++] = vInitModel[j].aa;
				break;
			case 5:
				vModel[indexModel++] = vInitModel[j].gm;
				break;					
			case 6:
				vModel[indexModel++] = vInitModel[j].az;
				break;					
			case 7:
				vModel[indexModel++] = vInitModel[j].S0;
				break;					
			case 8:
				vModel[indexModel++] = vInitModel[j].S1;
				break;					
			case 9:
				vModel[indexModel++] = vInitModel[j].mac;
				break;					
			case 10:
				vModel[indexModel++] = vInitModel[j].alfa;
				break;
			case 11: // NUMBER OF ITERATIONS
				vModel[indexModel++] = vNumIterPixel[j];
				break;
			case 12: // CHISQR 
				vModel[indexModel++] = vChisqrf[j];
				break;										
			default:
				break;			
			}
		}
	}

	long fpixel=1;
	if(displs>1){
		fpixel = (displs * naxes[2]) +1;
	}
	if ( fits_write_img(fptr, TFLOAT, displs, indexModel, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}
	/*fits_write_pix(fptr,TFLOAT, fpixel, indexModel,vModel,&status);
	if(status){
		fits_report_error(stderr, status);
		free(vModel);
		return 0;	
	}*/	
	free(vModel);
	if ( fits_close_file(fptr, &status) ){        
		printerror( status );
		return 0;
	}
	
	return 1;

}

/**
 * 
 * fixed = array with positions to write in the file, Positions are in the following order: 
 * [Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
 * */
/*int writeFitsImageModels(const char * fitsFile, int numRows, int numCols, Init_Model * vInitModel, PRECISION * vChisqrf, int * fixed, int addChiqr){

	fitsfile *fptr;       
   int status;
	int i, j, h; // indexes for loops
   long  fpixel, exposure;
	int indexModel = 0, sizeToCheck = NUMBER_PARAM_MODELS;


	int bitpix =  DOUBLE_IMG; 
   long naxis =   3;  
	long naxes[3] = { numRows, numCols, 0 };  

	for( i=0;i<NUMBER_PARAM_MODELS;i++){
		if(fixed[i]) naxes[2]++;
	}
	if(addChiqr){
		naxes[2]++;
		sizeToCheck++;
	}
   

   remove(fitsFile);               
   status = 0;         
   if (fits_create_file(&fptr, fitsFile, &status)) 
   	printerror( status );           
	

	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
		printerror( status );
		return 0;
	}

	PRECISION * vModel = calloc(naxes[0] * naxes[1] * naxes[2], sizeof(PRECISION));

	for( i=0;i<sizeToCheck;i++){
		if(i<NUMBER_PARAM_MODELS){
			if(fixed[i]){
				for( j=0;j<naxes[0];j++){
					for( h=0; h<naxes[1];h++){
						//[Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
						switch (i)
						{
						case 0:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].eta0;
							break;
						case 1:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].B;
							break;
						case 2:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].vlos;
							break;
						case 3:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].dopp;
							break;
						case 4:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].aa;
							break;
						case 5:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].gm;
							break;					
						case 6:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].az;
							break;					
						case 7:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].S0;
							break;					
						case 8:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].S1;
							break;					
						case 9:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].mac;
							break;					
						case 10:
							vModel[indexModel++] = vInitModel[( j*naxes[1]) + h].alfa;
							break;
						default:
							break;
						}
					}
				}
			}
		}
		else{  // add chisqr 
			for( j=0;j<naxes[0];j++){
				for( h=0; h<naxes[1];h++){		
					vModel[indexModel++] = vChisqrf[( j*naxes[1]) + h];
				}
			}
		}

	}

   fpixel = 1;                               
   //nelements = naxes[0] * naxes[1] * naxes[2];          


   if ( fits_write_img(fptr, TDOUBLE, fpixel, indexModel, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}

	// CLEAN MEMORY 
	free(vModel);


   exposure = 1500;
	if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
		"Total Exposure Time", &status) ){
		printerror( status );           
		return 0;
	}
	
	if ( fits_close_file(fptr, &status) ){        
		printerror( status );
		return 0;
	}
	
	return 1;

}

*/

int writeFitsImageProfiles(const char * fitsProfileFile, const char * fitsFileOrigin, FitsImage * image){

	fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
	int status = 0, ii = 1;
	int i, j, k, h; // indexes for loops
	int bitpix,naxis = 0, nkeys;
	long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	char card[FLEN_CARD];
	char keyname [FLEN_CARD];
	char value [FLEN_CARD];
	remove(fitsProfileFile);               /* Delete old file if it already exists */
	/* Open the input file and create output file */
   fits_open_file(&infptr, fitsFileOrigin, READONLY, &status);
   fits_create_file(&outfptr, fitsProfileFile, &status);
	if (status != 0) {    
		fits_report_error(stderr, status);
		return(status);
	}

	// read a maximun of 4 dimensions 
	fits_get_img_param(infptr, 4, &bitpix, &naxis, naxes, &status);
	// CREATE A NEW IMAGE 
   fits_create_img(outfptr, bitpix, naxis, naxes, &status);
	if (status) {
		fits_report_error(stderr, status);
		return(status);
	}
	/* copy all the user keywords (not the structural keywords) */
	fits_get_hdrspace(infptr, &nkeys, NULL, &status); 

	for (ii = 1; ii <= nkeys; ii++) {
		fits_read_record(infptr, ii, card, &status);
		fits_read_keyn(infptr,ii,keyname, value, NULL, &status);
		fits_update_card(outfptr, keyname,card, &status);
	}

	// CLOSE THE ORIGIN FILE, HE HAVE ALREADY THE INFORMATION OF KEYWORDS. 

	fits_close_file(infptr, &status);
	if (status){
		fits_report_error(stderr, status);
		return 0;
	}

	// ALLOCATE MEMORY TO WRITE THE IMAGE
	int numElemWrite = naxes[3]*naxes[2]*naxes[1]*naxes[0];

	float * outputImage = calloc(numElemWrite, sizeof(float));
	int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0;
	 
	int pos_lambda = image->pos_lambda;
	int pos_col = image->pos_col;
	int pos_row = image->pos_row;
	int pos_stokes_parameters = image->pos_stokes_parameters;
	
		
	for( i=0; i <naxes[3]; i++){
		for( j=0;j <naxes[2]; j++){
			for( k=0; k<naxes[1]; k++){
				for( h=0; h<naxes[0]; h++){
//	for( fpixel[3] = 1; fpixel[3]<=naxes[3];fpixel[3]++){
//		for( fpixel[2] = 1; fpixel[2]<=naxes[2];fpixel[2]++){
//			for( fpixel[1] = 1; fpixel[1]<=naxes[1];fpixel[1]++){
//				for( fpixel[0] = 1; fpixel[0]<=naxes[0]; fpixel[0]++){	
					// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
					switch (pos_lambda)
					{
						case 0:
							currentLambda = h;
							//currentLambda = fpixel[0]-1;
							break;
						case 1:
							currentLambda = k;
							//currentLambda = fpixel[1]-1;
							break;
						case 2:
							currentLambda = j;
							//currentLambda = fpixel[2]-1;
							break;
						case 3:
							currentLambda = i;
							//currentLambda = fpixel[3]-1;
							break;																						
					}
					switch (pos_stokes_parameters)
					{
						case 0:
							currentStokeParameter = h;
							//currentStokeParameter = fpixel[0]-1;
							break;
						case 1:
							currentStokeParameter = k;
							//currentStokeParameter = fpixel[1]-1;
							break;
						case 2:
							currentStokeParameter = j;
							//currentStokeParameter = fpixel[2]-1;
							break;
						case 3:
							currentStokeParameter = i;
							//currentStokeParameter = fpixel[3]-1;
							break;																						
					}
					switch (pos_row)
					{
						case 0:
							currentRow = h;
							//currentRow = fpixel[0]-1;
							break;
						case 1:
							currentRow = k;
							//currentRow = fpixel[1]-1;
							break;
						case 2:
							currentRow = j;
							//currentRow = fpixel[2]-1;
							break;
						case 3:
							currentRow = i;
							//currentRow = fpixel[3]-1;
							break;																						
					}
					switch (pos_col)
					{
						case 0:
							currentCol = h;
							//currentCol = fpixel[0]-1;
							break;
						case 1:
							currentCol = k;
							//currentCol = fpixel[1]-1;;
							break;
						case 2:
							currentCol = j;
							//currentCol = fpixel[2]-1;
							break;
						case 3:
							currentCol = i;
							//currentCol = fpixel[3]-1;
							break;																						
					}
					//double pixel = image->pixels[(currentRow*image->cols) + currentCol].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
					//fits_write_pix(outfptr, datatype, fpixel, 1, &pixel, &status);			
					outputImage[(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h] = image->pixels[(currentCol*image->rows) + currentRow].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
					//outputImage[numiter++] = image->pixels[(currentCol*image->rows) + currentRow].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
				}
			}
		}
	}

    /* write the array of unsigned integers to the FITS file */
	if ( fits_write_img(outfptr, TFLOAT, 1, numElemWrite, outputImage, &status) ){
		printerror( status );
		free(outputImage);
		return 0;
	}

	free(outputImage);
	fits_close_file(outfptr,  &status);
	if(status){
		printerror( status );
		return 0;
	}
	return 1;
}

int writeFitsImageModelsWithArray(char * fitsFile, int numRows, int numCols, PRECISION * eta0, PRECISION * B, PRECISION * vlos, PRECISION * dopp, PRECISION * aa, PRECISION * gm, PRECISION * az, PRECISION * S0, PRECISION * S1, PRECISION * mac, PRECISION * alfa, PRECISION * vChisqrf){

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
   int status;
	int i, j, h; // indexes for loops
   long  fpixel, nelements, exposure;

	int bitpix =  FLOAT_IMG; /* 16-bit unsigned short pixel values       */
   long naxis =   3;  /* 2-dimensional image                            */    
   long naxes[3] = { numRows, numCols, NUMBER_PARAM_MODELS };   /* Image of numRows X numCols x 10 parameters of model and chisqrf */

   remove(fitsFile);               /* Delete old file if it already exists */
   status = 0;         /* initialize status before calling fitsio routines */
   if (fits_create_file(&fptr, fitsFile, &status)) /* create new FITS file */
   	printerror( status );           /* call printerror if error occurs */
	
	 /* write the required keywords for the primary array image.     */
    /* Since bitpix = FLOAT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = -32 (float) .Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
		printerror( status );
		return 0;
	}

	float * vModel = malloc(numRows * numCols * NUMBER_PARAM_MODELS);

	int indexModel = 0;
	for( i=0;i<NUMBER_PARAM_MODELS;i++){
		for( j=0;j<numCols;j++){
			for( h=0; h<numRows;h++){
				switch (i)
				{
				case 0:
					vModel[indexModel++] = B[( j*numRows) + h];
					break;
				case 1:
					vModel[indexModel++] = gm[( j*numRows) + h];
					break;
				case 2:
					vModel[indexModel++] = az[( j*numRows) + h];
					break;
				case 3:
					vModel[indexModel++] = eta0[( j*numRows) + h];
					break;
				case 4:
					vModel[indexModel++] = dopp[( j*numRows) + h];
					break;
				case 5:
					vModel[indexModel++] = aa[( j*numRows) + h];
					break;					
				case 6:
					vModel[indexModel++] = vlos[( j*numRows) + h];
					break;					
				case 7:
					vModel[indexModel++] = alfa[( j*numRows) + h];
					break;					
				case 8:
					vModel[indexModel++] = S0[( j*numRows) + h];
					break;					
				case 9:
					vModel[indexModel++] = S1[( j*numRows) + h];
					break;					
				case 10:
					vModel[indexModel++] = mac[( j*numRows) + h];
					break;										
				case 11: // READ FROM CHISQR
					vModel[indexModel++] = vChisqrf[( j*numRows) + h];
					break;					
				default:
					break;
				}
			}
		}
	}

   fpixel = 1;                               /* first pixel to write      */
   nelements = naxes[0] * naxes[1] * naxes[2];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */
   if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}

	// CLEAN MEMORY 
	free(vModel);

	    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    exposure = 1500;
	if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
		"Total Exposure Time", &status) ){
		printerror( status );           
		return 0;
	}

	if ( fits_close_file(fptr, &status) ){                /* close the file */
		printerror( status );
		return 0;
	}
	
	return 1;

}


int readSizeImageSpectro(const char * fitsFile, int * numRows, int * numCols){
	fitsfile *fptr;   
   int status = 0;   
	FitsImage * image =  malloc(sizeof(FitsImage));
	PRECISION nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, naxis, anynul, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; 
	char comment[FLEN_CARD];

	if (!fits_open_file(&fptr, fitsFile, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);
		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;
			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_lambda; 
			int pos_row;
			int pos_col;
			int pos_stokes_parameters;
			// LAMBDA POSITION
			if(strcmp(image->ctype_1,CTYPE_WAVE)==0) pos_lambda = 0;
			if(strcmp(image->ctype_2,CTYPE_WAVE)==0) pos_lambda = 1;
			if(strcmp(image->ctype_3,CTYPE_WAVE)==0) pos_lambda = 2;
			if(strcmp(image->ctype_4,CTYPE_WAVE)==0) pos_lambda = 3;

			// HPLN TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)==0) pos_row = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLN_TAN)==0) pos_row = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLN_TAN)==0) pos_row = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLN_TAN)==0) pos_row = 3;

			// HPLT TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLT_TAN)==0) pos_col = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)==0) pos_col = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLT_TAN)==0) pos_col = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLT_TAN)==0) pos_col = 3;			

			// Stokes paramter position , 
			if(strcmp(image->ctype_1,CTYPE_STOKES)==0) pos_stokes_parameters = 0;
			if(strcmp(image->ctype_2,CTYPE_STOKES)==0) pos_stokes_parameters = 1;
			if(strcmp(image->ctype_3,CTYPE_STOKES)==0) pos_stokes_parameters = 2;
			if(strcmp(image->ctype_4,CTYPE_STOKES)==0) pos_stokes_parameters = 3;

			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
				*numRows=naxes[pos_row];
				*numCols=naxes[pos_col];
				free(image);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return 0;
				}
			}
			else{
				printf("\n************ ERROR getting the image from the fits file:  %s",fitsFile);
				fits_close_file(fptr, &status);
				free(image);
				if (status){
					fits_report_error(stderr, status);
					return 0;
				}			
				return 0;				
			}

		}
		else
		{
			printf("\n************ ERROR: Fits file: %s could not be a fits image.",fitsFile);
			fits_close_file(fptr, &status);
			free(image);
			if (status){
				fits_report_error(stderr, status);
				return 0;
			}			
			return 0;
		}
		
	}
	else{
		printf("\n************ ERROR openning fits file %s",fitsFile);
		fits_close_file(fptr, &status);
		free(image);
		if (status){
			fits_report_error(stderr, status);
			return 0;
		}			
	}
	free(image);
	return 1;

}

/*--------------------------------------------------------------------------*/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       //exit( status );    /* terminate the program, returning error status */
    }
    return;
}



