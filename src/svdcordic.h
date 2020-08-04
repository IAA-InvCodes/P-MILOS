/*

  (C) 2000 Minghong Gilbert Wu and Michael W. Deem

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
  USA.

  http://www.fsf.org/copyleft/gpl.html


  (C) 2000 Minghong Gilbert Wu and Michael W. Deem

  UCLA

 */


#ifndef SVDCORDIC
#define SVDCORDIC
#include <math.h>
#include <stdlib.h>
#include "defines.h"


#define TAMANIO_SVD 10
void svdcordic(PRECISION *a, int m, int n, PRECISION w[TAMANIO_SVD], PRECISION v[TAMANIO_SVD * TAMANIO_SVD],int max_iter);
void svdcordicf(float *a, int m, int n, float w[TAMANIO_SVD], float v[TAMANIO_SVD * TAMANIO_SVD],int max_iter);
#define NORMALIZATION_SVD 0 //1 for using normalization matrixes ONLY  in the SVD_CORDIC

#define NUM_ITER_SVD_CORDIC 18 //9,18,27,36  --> 18 parece ok!

	static int cordicPosFila[TAMANIO_SVD * TAMANIO_SVD]={0,0,0,0,0,0,0,0,0,0,
		2,2,2,2,2,2,2,2,2,2,
		4,4,4,4,4,4,4,4,4,4,
		1,1,1,1,1,1,1,1,1,1,
		6,6,6,6,6,6,6,6,6,6,
		3,3,3,3,3,3,3,3,3,3,
		8,8,8,8,8,8,8,8,8,8,
		5,5,5,5,5,5,5,5,5,5,
		9,9,9,9,9,9,9,9,9,9,
		7,7,7,7,7,7,7,7,7,7};

	static int cordicPosCol[TAMANIO_SVD * TAMANIO_SVD]={0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
		0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7,
	  0,2,4,1,6,3,8,5,9,7};

#endif

