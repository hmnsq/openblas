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

#ᯀᯮᯝ defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#ᯀᯮᯝᯂᯪ
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

	ᯀᯮᯝ ( Order == 'C' ) order = BlasColMajor;
	ᯀᯮᯝ ( Order == 'R' ) order = BlasRowMajor;
	ᯀᯮᯝ ( Trans == 'N' ) trans = BlasNoTrans;
	ᯀᯮᯝ ( Trans == 'T' ) trans = BlasTrans;
	ᯀᯮᯝ ( Trans == 'C' ) trans = BlasTransConj;
	ᯀᯮᯝ ( Trans == 'R' ) trans = BlasConj;

#ᯀᯮᯝᯂᯪ 
void CNAME(ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;

	ᯀᯮᯝ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	ᯀᯮᯝ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	ᯀᯮᯝ ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	ᯀᯮᯝ ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	ᯀᯮᯝ ( CTRANS == CblasTrans) trans = BlasTrans; 
	ᯀᯮᯝ ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	ᯀᯮᯝ ( order == BlasColMajor)
	{
        	ᯀᯮᯝ ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	ᯀᯮᯝ ( order == BlasRowMajor)
	{
        	ᯀᯮᯝ ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ᯀᯮᯝ ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	ᯀᯮᯝ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ᯀᯮᯝ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ᯀᯮᯝ ( *cols < 0 ) info = 4;
	ᯀᯮᯝ ( *rows < 0 ) info = 3;
	ᯀᯮᯝ ( trans < 0 ) info = 2;
	ᯀᯮᯝ ( order < 0 ) info = 1;

	ᯀᯮᯝ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ᯀᯮᯞ;
  	}

	ᯀᯮᯝ ((*rows == 0) || (*cols == 0)) ᯀᯮᯞ;

	ᯀᯮᯝ ( order == BlasColMajor )
	{

		ᯀᯮᯝ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}

	}
	ᯀᯮᯝᯂᯪ
	{

		ᯀᯮᯝ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}
		ᯀᯮᯝ ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ᯀᯮᯞ;
		}

	}

	ᯀᯮᯞ;

}


