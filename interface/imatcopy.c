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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ଯଦି defined(DOUBLE)
#define ERROR_NAME "DIMATCOPY"
#ଅନ୍ୟଥା
#define ERROR_NAME "SIMATCOPY"
#endif

#define BlasRowMajor 0
#define BlasColMajor 1
#define BlasNoTrans  0
#define BlasTrans    1

#undef malloc
#undef ମୁକ୍ତ

/* Enables the New IMATCOPY code with inplace operation if lda == ldb   */
#define NEW_IMATCOPY

OPENBLAS_EXPORT

#ifndef CBLAS
ଶୂନ୍ୟ NAME( ଅକ୍ଷର* ORDER, ଅକ୍ଷର* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, blasint *ldb)
{

	ଅକ୍ଷର Order, Trans;
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	ଯଦି ( Order == 'C' ) order = BlasColMajor;
	ଯଦି ( Order == 'R' ) order = BlasRowMajor;
	ଯଦି ( Trans == 'N' ) trans = BlasNoTrans;
	ଯଦି ( Trans == 'R' ) trans = BlasNoTrans;
	ଯଦି ( Trans == 'T' ) trans = BlasTrans;
	ଯଦି ( Trans == 'C' ) trans = BlasTrans;
#ଅନ୍ୟଥା 
ଶୂନ୍ୟ CNAME( ଗଣନା CBLAS_ORDER CORDER, ଗଣନା CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, blasint cldb)
{
	ପୂର୍ଣ୍ଣସଂଖ୍ୟା order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;
	blasint *lda, *ldb, *rows, *cols; 
	FLOAT *alpha; 

	ଯଦି ( CORDER == CblasColMajor) order = BlasColMajor; 
	ଯଦି ( CORDER == CblasRowMajor) order = BlasRowMajor; 
	ଯଦି ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans) trans = BlasNoTrans; 
	ଯଦି ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans  ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	alpha = &calpha; 
	lda = &clda; 
	ldb = &cldb; 	
#endif 

	ଯଦି ( order == BlasColMajor)
	{
        	ଯଦି ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 8;
        	ଯଦି ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 8;
	}
	ଯଦି ( order == BlasRowMajor)
	{
        	ଯଦି ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 8;
        	ଯଦି ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 8;
	}

	ଯଦି ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ଯଦି ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ଯଦି ( *cols < 0 ) info = 4;
	ଯଦି ( *rows < 0 ) info = 3;
	ଯଦି ( trans < 0 ) info = 2;
	ଯଦି ( order < 0 ) info = 1;

	ଯଦି (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    		ଫେରନ୍ତୁ;
  	}

	ଯଦି ((*rows == 0) || (*cols == 0)) ଫେରନ୍ତୁ;

#ifdef NEW_IMATCOPY
    ଯଦି ( *lda == *ldb ) {
        ଯଦି ( order == BlasColMajor )
        {
            ଯଦି ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda );
                ଫେରନ୍ତୁ;
            }
            ଅନ୍ୟଥା ଯଦି ( *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda );
                ଫେରନ୍ତୁ;
            }
        }
        ଅନ୍ୟଥା
        {
            ଯଦି ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda );
                ଫେରନ୍ତୁ;
            }
            ଅନ୍ୟଥା ଯଦି ( *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda );
                ଫେରନ୍ତୁ;
            }
        }
    }
#endif

	ଯଦି ( *rows >  *cols )
            msize = (size_t)(*rows) * (*ldb)  * ଆକାର(FLOAT);
    ଅନ୍ୟଥା
            msize = (size_t)(*cols) * (*ldb)  * ଆକାର(FLOAT);

	b = malloc(msize);
	ଯଦି ( b == ଖାଲି )
	{
		printf("Memory alloc failed in imatcopy\n");
		exit(1);
	}

	ଯଦି ( order == BlasColMajor )
	{
		ଯଦି ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0 , b, *rows, a, *ldb );
		}
		ଅନ୍ୟଥା
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
	}
	ଅନ୍ୟଥା
	{
		ଯଦି ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *cols );
			OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, b, *cols, a, *ldb );
		}
		ଅନ୍ୟଥା
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *rows );
			OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, b, *rows, a, *ldb );
		}
	}

	ମୁକ୍ତ(b);
	ଫେରନ୍ତୁ;

}


