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
 * 2014/06/09 Saar
***********************************************************/

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#اگر defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#وگرنه
#define ERROR_NAME "COMATCOPY"
#endif

#define BlasRowMajor     0
#define BlasColMajor     1
#define BlasNoTrans      0
#define BlasTrans        1
#define BlasTransConj    2
#define BlasConj         3

OPENBLAS_EXPORT

#ifndef CBLAS
تهی NAME( نویسه* ORDER, نویسه* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	نویسه Order, Trans;
	صحیح order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	اگر ( Order == 'C' ) order = BlasColMajor;
	اگر ( Order == 'R' ) order = BlasRowMajor;
	اگر ( Trans == 'N' ) trans = BlasNoTrans;
	اگر ( Trans == 'T' ) trans = BlasTrans;
	اگر ( Trans == 'C' ) trans = BlasTransConj;
	اگر ( Trans == 'R' ) trans = BlasConj;

#وگرنه 
تهی CNAME(شمارش CBLAS_ORDER CORDER, شمارش CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	صحیح order=-1,trans=-1;
	blasint info = -1;

	اگر ( CORDER == CblasColMajor ) order = BlasColMajor; 
	اگر ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	اگر ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	اگر ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	اگر ( CTRANS == CblasTrans) trans = BlasTrans; 
	اگر ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	اگر ( order == BlasColMajor)
	{
        	اگر ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	اگر ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	اگر ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	اگر ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	اگر ( order == BlasRowMajor)
	{
        	اگر ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	اگر ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	اگر ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	اگر ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	اگر ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	اگر ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	اگر ( *cols < 0 ) info = 4;
	اگر ( *rows < 0 ) info = 3;
	اگر ( trans < 0 ) info = 2;
	اگر ( order < 0 ) info = 1;

	اگر (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    		بازگشت;
  	}

	اگر ((*rows == 0) || (*cols == 0)) بازگشت;

	اگر ( order == BlasColMajor )
	{

		اگر ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}

	}
	وگرنه
	{

		اگر ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}
		اگر ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			بازگشت;
		}

	}

	بازگشت;

}


