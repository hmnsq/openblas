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

#ᬬᬾᬦ᭄ defined(DOUBLE)
#define ERROR_NAME "ZIMATCOPY"
#ᬮᬶᬬ
#define ERROR_NAME "CIMATCOPY"
#endif

#define BlasRowMajor     0
#define BlasColMajor     1
#define BlasNoTrans      0
#define BlasTrans        1
#define BlasTransConj    2
#define BlasConj         3

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

	ᬬᬾᬦ᭄ ( Order == 'C' ) order = BlasColMajor;
	ᬬᬾᬦ᭄ ( Order == 'R' ) order = BlasRowMajor;
	ᬬᬾᬦ᭄ ( Trans == 'N' ) trans = BlasNoTrans;
	ᬬᬾᬦ᭄ ( Trans == 'T' ) trans = BlasTrans;
	ᬬᬾᬦ᭄ ( Trans == 'C' ) trans = BlasTransConj;
	ᬬᬾᬦ᭄ ( Trans == 'R' ) trans = BlasConj;

#ᬮᬶᬬ 
void CNAME( ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT *alpha, FLOAT *a, blasint clda, blasint cldb)
{

	blasint *rows, *cols, *lda, *ldb; 
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	ᬬᬾᬦ᭄ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	ᬬᬾᬦ᭄ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	ᬬᬾᬦ᭄ ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	ᬬᬾᬦ᭄ ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	ᬬᬾᬦ᭄ ( CTRANS == CblasTrans) trans = BlasTrans; 
	ᬬᬾᬦ᭄ ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif

	ᬬᬾᬦ᭄ ( order == BlasColMajor)
	{
        	ᬬᬾᬦ᭄ ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	ᬬᬾᬦ᭄ ( order == BlasRowMajor)
	{
        	ᬬᬾᬦ᭄ ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᬬᬾᬦ᭄ ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	ᬬᬾᬦ᭄ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ᬬᬾᬦ᭄ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ᬬᬾᬦ᭄ ( *cols < 0 ) info = 4;
	ᬬᬾᬦ᭄ ( *rows < 0 ) info = 3;
	ᬬᬾᬦ᭄ ( trans < 0 ) info = 2;
	ᬬᬾᬦ᭄ ( order < 0 ) info = 1;

	ᬬᬾᬦ᭄ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ᬯᬮᬶ;
  	}

	ᬬᬾᬦ᭄ ((*rows == 0) || (*cols == 0)) ᬯᬮᬶ;

#ifdef NEW_IMATCOPY
    ᬬᬾᬦ᭄ (*lda == *ldb ) {
        ᬬᬾᬦ᭄ ( order == BlasColMajor )
        {

            ᬬᬾᬦ᭄ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasConj )
            {
                IMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasTrans && *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasTransConj && *rows == *cols )
            {
                IMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }

        }
        ᬮᬶᬬ
        {

            ᬬᬾᬦ᭄ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasConj )
            {
                IMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasTrans && *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }
            ᬬᬾᬦ᭄ ( trans == BlasTransConj && *rows == *cols )
            {
                IMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ᬯᬮᬶ;
            }

        }
    }
#endif

		ᬬᬾᬦ᭄ ( *rows >  *cols )
                msize = (size_t)(*rows) * (*ldb)  * parimaa_na(FLOAT) * 2;
        ᬮᬶᬬ
                msize = (size_t)(*cols) * (*ldb)  * parimaa_na(FLOAT) * 2;

	b = malloc(msize);
	ᬬᬾᬦ᭄ ( b == NULL )
	{
		printf("Memory alloc failed in zimatcopy\n");
		exit(1);
	}

	ᬬᬾᬦ᭄ ( order == BlasColMajor )
	{

		ᬬᬾᬦ᭄ ( trans == BlasNoTrans )
		{
	  		OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}

	}
	ᬮᬶᬬ
	{

		ᬬᬾᬦ᭄ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ᬮᬶᬬ ᬬᬾᬦ᭄ ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}

	}

	free(b);
	ᬯᬮᬶ;

}


