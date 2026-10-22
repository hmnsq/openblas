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

#ᰀᰧ defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#ᰀᰧᰀ
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
void NAME( var_na* ORDER, var_na* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	var_na Order, Trans;
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	ᰀᰧ ( Order == 'C' ) order = BlasColMajor;
	ᰀᰧ ( Order == 'R' ) order = BlasRowMajor;
	ᰀᰧ ( Trans == 'N' ) trans = BlasNoTrans;
	ᰀᰧ ( Trans == 'T' ) trans = BlasTrans;
	ᰀᰧ ( Trans == 'C' ) trans = BlasTransConj;
	ᰀᰧ ( Trans == 'R' ) trans = BlasConj;

#ᰀᰧᰀ 
void CNAME(ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;

	ᰀᰧ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	ᰀᰧ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	ᰀᰧ ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	ᰀᰧ ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	ᰀᰧ ( CTRANS == CblasTrans) trans = BlasTrans; 
	ᰀᰧ ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	ᰀᰧ ( order == BlasColMajor)
	{
        	ᰀᰧ ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᰀᰧ ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᰀᰧ ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᰀᰧ ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	ᰀᰧ ( order == BlasRowMajor)
	{
        	ᰀᰧ ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᰀᰧ ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᰀᰧ ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᰀᰧ ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	ᰀᰧ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ᰀᰧ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ᰀᰧ ( *cols < 0 ) info = 4;
	ᰀᰧ ( *rows < 0 ) info = 3;
	ᰀᰧ ( trans < 0 ) info = 2;
	ᰀᰧ ( order < 0 ) info = 1;

	ᰀᰧ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ᰀᰧᰌ;
  	}

	ᰀᰧ ((*rows == 0) || (*cols == 0)) ᰀᰧᰌ;

	ᰀᰧ ( order == BlasColMajor )
	{

		ᰀᰧ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}

	}
	ᰀᰧᰀ
	{

		ᰀᰧ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}
		ᰀᰧ ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᰀᰧᰌ;
		}

	}

	ᰀᰧᰌ;

}


