/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/***********************************************************
 * 2014-06-10 Saar
 * 2015-09-07 grisuthedragon 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ถ้า defined(DOUBLE)
#define ERROR_NAME "DIMATCOPY"
#มิฉะนั้น
#define ERROR_NAME "SIMATCOPY"
#endif

#define BlasRowMajor 0
#define BlasColMajor 1
#define BlasNoTrans  0
#define BlasTrans    1

#undef malloc
#undef free

/* Enables the New IMATCOPY code with inplace operation if lda == ldb   */
#define NEW_IMATCOPY

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME( var_na* ORDER, var_na* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, blasint *ldb)
{

	var_na Order, Trans;
	จำนวนเต็ม order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	ถ้า ( Order == 'C' ) order = BlasColMajor;
	ถ้า ( Order == 'R' ) order = BlasRowMajor;
	ถ้า ( Trans == 'N' ) trans = BlasNoTrans;
	ถ้า ( Trans == 'R' ) trans = BlasNoTrans;
	ถ้า ( Trans == 'T' ) trans = BlasTrans;
	ถ้า ( Trans == 'C' ) trans = BlasTrans;
#มิฉะนั้น 
void CNAME( ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, blasint cldb)
{
	จำนวนเต็ม order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;
	blasint *lda, *ldb, *rows, *cols; 
	FLOAT *alpha; 

	ถ้า ( CORDER == CblasColMajor) order = BlasColMajor; 
	ถ้า ( CORDER == CblasRowMajor) order = BlasRowMajor; 
	ถ้า ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans) trans = BlasNoTrans; 
	ถ้า ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans  ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	alpha = &calpha; 
	lda = &clda; 
	ldb = &cldb; 	
#endif 

	ถ้า ( order == BlasColMajor)
	{
        	ถ้า ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 8;
        	ถ้า ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 8;
	}
	ถ้า ( order == BlasRowMajor)
	{
        	ถ้า ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 8;
        	ถ้า ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 8;
	}

	ถ้า ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ถ้า ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ถ้า ( *cols < 0 ) info = 4;
	ถ้า ( *rows < 0 ) info = 3;
	ถ้า ( trans < 0 ) info = 2;
	ถ้า ( order < 0 ) info = 1;

	ถ้า (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		คืน;
  	}

	ถ้า ((*rows == 0) || (*cols == 0)) คืน;

#ifdef NEW_IMATCOPY
    ถ้า ( *lda == *ldb ) {
        ถ้า ( order == BlasColMajor )
        {
            ถ้า ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda );
                คืน;
            }
            มิฉะนั้น ถ้า ( *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda );
                คืน;
            }
        }
        มิฉะนั้น
        {
            ถ้า ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda );
                คืน;
            }
            มิฉะนั้น ถ้า ( *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda );
                คืน;
            }
        }
    }
#endif

	ถ้า ( *rows >  *cols )
            msize = (size_t)(*rows) * (*ldb)  * parimaa_na(FLOAT);
    มิฉะนั้น
            msize = (size_t)(*cols) * (*ldb)  * parimaa_na(FLOAT);

	b = malloc(msize);
	ถ้า ( b == NULL )
	{
		printf("Memory alloc failed in imatcopy\n");
		exit(1);
	}

	ถ้า ( order == BlasColMajor )
	{
		ถ้า ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0 , b, *rows, a, *ldb );
		}
		มิฉะนั้น
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
	}
	มิฉะนั้น
	{
		ถ้า ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
		มิฉะนั้น
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, b, *rows, a, *ldb );
		}
	}

	free(b);
	คืน;

}


