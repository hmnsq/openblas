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

#nếu defined(DOUBLE)
#define ERROR_NAME "DIMATCOPY"
#không_thì
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
void NAME( char* ORDER, char* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, blasint *ldb)
{

	char Order, Trans;
	số_nguyên order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	nếu ( Order == 'C' ) order = BlasColMajor;
	nếu ( Order == 'R' ) order = BlasRowMajor;
	nếu ( Trans == 'N' ) trans = BlasNoTrans;
	nếu ( Trans == 'R' ) trans = BlasNoTrans;
	nếu ( Trans == 'T' ) trans = BlasTrans;
	nếu ( Trans == 'C' ) trans = BlasTrans;
#không_thì 
void CNAME( enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, blasint cldb)
{
	số_nguyên order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;
	blasint *lda, *ldb, *rows, *cols; 
	FLOAT *alpha; 

	nếu ( CORDER == CblasColMajor) order = BlasColMajor; 
	nếu ( CORDER == CblasRowMajor) order = BlasRowMajor; 
	nếu ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans) trans = BlasNoTrans; 
	nếu ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans  ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	alpha = &calpha; 
	lda = &clda; 
	ldb = &cldb; 	
#endif 

	nếu ( order == BlasColMajor)
	{
        	nếu ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 8;
        	nếu ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 8;
	}
	nếu ( order == BlasRowMajor)
	{
        	nếu ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 8;
        	nếu ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 8;
	}

	nếu ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	nếu ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	nếu ( *cols < 0 ) info = 4;
	nếu ( *rows < 0 ) info = 3;
	nếu ( trans < 0 ) info = 2;
	nếu ( order < 0 ) info = 1;

	nếu (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    		trả_về;
  	}

	nếu ((*rows == 0) || (*cols == 0)) trả_về;

#ifdef NEW_IMATCOPY
    nếu ( *lda == *ldb ) {
        nếu ( order == BlasColMajor )
        {
            nếu ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda );
                trả_về;
            }
            không_thì nếu ( *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda );
                trả_về;
            }
        }
        không_thì
        {
            nếu ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda );
                trả_về;
            }
            không_thì nếu ( *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda );
                trả_về;
            }
        }
    }
#endif

	nếu ( *rows >  *cols )
            msize = (size_t)(*rows) * (*ldb)  * sizeof(FLOAT);
    không_thì
            msize = (size_t)(*cols) * (*ldb)  * sizeof(FLOAT);

	b = malloc(msize);
	nếu ( b == NULL )
	{
		printf("Memory alloc failed in imatcopy\n");
		exit(1);
	}

	nếu ( order == BlasColMajor )
	{
		nếu ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0 , b, *rows, a, *ldb );
		}
		không_thì
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
	}
	không_thì
	{
		nếu ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
		không_thì
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, b, *rows, a, *ldb );
		}
	}

	free(b);
	trả_về;

}


