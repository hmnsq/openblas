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

#ᱡᱩᱫᱤ defined(DOUBLE)
#define ERROR_NAME "DIMATCOPY"
#ᱵᱟᱝ
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
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	ᱡᱩᱫᱤ ( Order == 'C' ) order = BlasColMajor;
	ᱡᱩᱫᱤ ( Order == 'R' ) order = BlasRowMajor;
	ᱡᱩᱫᱤ ( Trans == 'N' ) trans = BlasNoTrans;
	ᱡᱩᱫᱤ ( Trans == 'R' ) trans = BlasNoTrans;
	ᱡᱩᱫᱤ ( Trans == 'T' ) trans = BlasTrans;
	ᱡᱩᱫᱤ ( Trans == 'C' ) trans = BlasTrans;
#ᱵᱟᱝ 
void CNAME( ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, blasint cldb)
{
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;
	blasint *lda, *ldb, *rows, *cols; 
	FLOAT *alpha; 

	ᱡᱩᱫᱤ ( CORDER == CblasColMajor) order = BlasColMajor; 
	ᱡᱩᱫᱤ ( CORDER == CblasRowMajor) order = BlasRowMajor; 
	ᱡᱩᱫᱤ ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans) trans = BlasNoTrans; 
	ᱡᱩᱫᱤ ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans  ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	alpha = &calpha; 
	lda = &clda; 
	ldb = &cldb; 	
#endif 

	ᱡᱩᱫᱤ ( order == BlasColMajor)
	{
        	ᱡᱩᱫᱤ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 8;
        	ᱡᱩᱫᱤ ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 8;
	}
	ᱡᱩᱫᱤ ( order == BlasRowMajor)
	{
        	ᱡᱩᱫᱤ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 8;
        	ᱡᱩᱫᱤ ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 8;
	}

	ᱡᱩᱫᱤ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ᱡᱩᱫᱤ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ᱡᱩᱫᱤ ( *cols < 0 ) info = 4;
	ᱡᱩᱫᱤ ( *rows < 0 ) info = 3;
	ᱡᱩᱫᱤ ( trans < 0 ) info = 2;
	ᱡᱩᱫᱤ ( order < 0 ) info = 1;

	ᱡᱩᱫᱤ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ᱨᱩᱣᱟᱹᱲ;
  	}

	ᱡᱩᱫᱤ ((*rows == 0) || (*cols == 0)) ᱨᱩᱣᱟᱹᱲ;

#ifdef NEW_IMATCOPY
    ᱡᱩᱫᱤ ( *lda == *ldb ) {
        ᱡᱩᱫᱤ ( order == BlasColMajor )
        {
            ᱡᱩᱫᱤ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda );
                ᱨᱩᱣᱟᱹᱲ;
            }
            ᱵᱟᱝ ᱡᱩᱫᱤ ( *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda );
                ᱨᱩᱣᱟᱹᱲ;
            }
        }
        ᱵᱟᱝ
        {
            ᱡᱩᱫᱤ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda );
                ᱨᱩᱣᱟᱹᱲ;
            }
            ᱵᱟᱝ ᱡᱩᱫᱤ ( *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda );
                ᱨᱩᱣᱟᱹᱲ;
            }
        }
    }
#endif

	ᱡᱩᱫᱤ ( *rows >  *cols )
            msize = (size_t)(*rows) * (*ldb)  * parimaa_na(FLOAT);
    ᱵᱟᱝ
            msize = (size_t)(*cols) * (*ldb)  * parimaa_na(FLOAT);

	b = malloc(msize);
	ᱡᱩᱫᱤ ( b == NULL )
	{
		printf("Memory alloc failed in imatcopy\n");
		exit(1);
	}

	ᱡᱩᱫᱤ ( order == BlasColMajor )
	{
		ᱡᱩᱫᱤ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0 , b, *rows, a, *ldb );
		}
		ᱵᱟᱝ
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
	}
	ᱵᱟᱝ
	{
		ᱡᱩᱫᱤ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
		ᱵᱟᱝ
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, b, *rows, a, *ldb );
		}
	}

	free(b);
	ᱨᱩᱣᱟᱹᱲ;

}


