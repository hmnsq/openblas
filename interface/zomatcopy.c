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

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#만약 defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#아니면
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
void NAME( char* ORDER, char* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	char Order, Trans;
	정수 order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	만약 ( Order == 'C' ) order = BlasColMajor;
	만약 ( Order == 'R' ) order = BlasRowMajor;
	만약 ( Trans == 'N' ) trans = BlasNoTrans;
	만약 ( Trans == 'T' ) trans = BlasTrans;
	만약 ( Trans == 'C' ) trans = BlasTransConj;
	만약 ( Trans == 'R' ) trans = BlasConj;

#아니면 
void CNAME(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	정수 order=-1,trans=-1;
	blasint info = -1;

	만약 ( CORDER == CblasColMajor ) order = BlasColMajor; 
	만약 ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	만약 ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	만약 ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	만약 ( CTRANS == CblasTrans) trans = BlasTrans; 
	만약 ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	만약 ( order == BlasColMajor)
	{
        	만약 ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	만약 ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	만약 ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	만약 ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	만약 ( order == BlasRowMajor)
	{
        	만약 ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	만약 ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	만약 ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	만약 ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	만약 ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	만약 ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	만약 ( *cols < 0 ) info = 4;
	만약 ( *rows < 0 ) info = 3;
	만약 ( trans < 0 ) info = 2;
	만약 ( order < 0 ) info = 1;

	만약 (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    		반환;
  	}

	만약 ((*rows == 0) || (*cols == 0)) 반환;

	만약 ( order == BlasColMajor )
	{

		만약 ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}

	}
	아니면
	{

		만약 ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}
		만약 ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			반환;
		}

	}

	반환;

}


