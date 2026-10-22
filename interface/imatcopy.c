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

#གལ་ཏེ defined(DOUBLE)
#define ERROR_NAME "DIMATCOPY"
#དེ་མིན
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
	ཧྲིལ order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	གལ་ཏེ ( Order == 'C' ) order = BlasColMajor;
	གལ་ཏེ ( Order == 'R' ) order = BlasRowMajor;
	གལ་ཏེ ( Trans == 'N' ) trans = BlasNoTrans;
	གལ་ཏེ ( Trans == 'R' ) trans = BlasNoTrans;
	གལ་ཏེ ( Trans == 'T' ) trans = BlasTrans;
	གལ་ཏེ ( Trans == 'C' ) trans = BlasTrans;
#དེ་མིན 
void CNAME( ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, blasint cldb)
{
	ཧྲིལ order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;
	blasint *lda, *ldb, *rows, *cols; 
	FLOAT *alpha; 

	གལ་ཏེ ( CORDER == CblasColMajor) order = BlasColMajor; 
	གལ་ཏེ ( CORDER == CblasRowMajor) order = BlasRowMajor; 
	གལ་ཏེ ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans) trans = BlasNoTrans; 
	གལ་ཏེ ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans  ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	alpha = &calpha; 
	lda = &clda; 
	ldb = &cldb; 	
#endif 

	གལ་ཏེ ( order == BlasColMajor)
	{
        	གལ་ཏེ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 8;
        	གལ་ཏེ ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 8;
	}
	གལ་ཏེ ( order == BlasRowMajor)
	{
        	གལ་ཏེ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 8;
        	གལ་ཏེ ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 8;
	}

	གལ་ཏེ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	གལ་ཏེ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	གལ་ཏེ ( *cols < 0 ) info = 4;
	གལ་ཏེ ( *rows < 0 ) info = 3;
	གལ་ཏེ ( trans < 0 ) info = 2;
	གལ་ཏེ ( order < 0 ) info = 1;

	གལ་ཏེ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ཕྱིར;
  	}

	གལ་ཏེ ((*rows == 0) || (*cols == 0)) ཕྱིར;

#ifdef NEW_IMATCOPY
    གལ་ཏེ ( *lda == *ldb ) {
        གལ་ཏེ ( order == BlasColMajor )
        {
            གལ་ཏེ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda );
                ཕྱིར;
            }
            དེ་མིན གལ་ཏེ ( *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda );
                ཕྱིར;
            }
        }
        དེ་མིན
        {
            གལ་ཏེ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda );
                ཕྱིར;
            }
            དེ་མིན གལ་ཏེ ( *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda );
                ཕྱིར;
            }
        }
    }
#endif

	གལ་ཏེ ( *rows >  *cols )
            msize = (size_t)(*rows) * (*ldb)  * parimaa_na(FLOAT);
    དེ་མིན
            msize = (size_t)(*cols) * (*ldb)  * parimaa_na(FLOAT);

	b = malloc(msize);
	གལ་ཏེ ( b == NULL )
	{
		printf("Memory alloc failed in imatcopy\n");
		exit(1);
	}

	གལ་ཏེ ( order == BlasColMajor )
	{
		གལ་ཏེ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0 , b, *rows, a, *ldb );
		}
		དེ་མིན
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
	}
	དེ་མིན
	{
		གལ་ཏེ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
		དེ་མིན
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, b, *rows, a, *ldb );
		}
	}

	free(b);
	ཕྱིར;

}


