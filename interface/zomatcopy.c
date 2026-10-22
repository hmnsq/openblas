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

#ꦪꦼꦤ꧀ defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#ꦭꦶꦪ
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

	ꦪꦼꦤ꧀ ( Order == 'C' ) order = BlasColMajor;
	ꦪꦼꦤ꧀ ( Order == 'R' ) order = BlasRowMajor;
	ꦪꦼꦤ꧀ ( Trans == 'N' ) trans = BlasNoTrans;
	ꦪꦼꦤ꧀ ( Trans == 'T' ) trans = BlasTrans;
	ꦪꦼꦤ꧀ ( Trans == 'C' ) trans = BlasTransConj;
	ꦪꦼꦤ꧀ ( Trans == 'R' ) trans = BlasConj;

#ꦭꦶꦪ 
void CNAME(ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	puur_naa_nka order=-1,trans=-1;
	blasint info = -1;

	ꦪꦼꦤ꧀ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	ꦪꦼꦤ꧀ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	ꦪꦼꦤ꧀ ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	ꦪꦼꦤ꧀ ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	ꦪꦼꦤ꧀ ( CTRANS == CblasTrans) trans = BlasTrans; 
	ꦪꦼꦤ꧀ ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	ꦪꦼꦤ꧀ ( order == BlasColMajor)
	{
        	ꦪꦼꦤ꧀ ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	ꦪꦼꦤ꧀ ( order == BlasRowMajor)
	{
        	ꦪꦼꦤ꧀ ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ꦪꦼꦤ꧀ ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	ꦪꦼꦤ꧀ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ꦪꦼꦤ꧀ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ꦪꦼꦤ꧀ ( *cols < 0 ) info = 4;
	ꦪꦼꦤ꧀ ( *rows < 0 ) info = 3;
	ꦪꦼꦤ꧀ ( trans < 0 ) info = 2;
	ꦪꦼꦤ꧀ ( order < 0 ) info = 1;

	ꦪꦼꦤ꧀ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ꦮꦭꦶ;
  	}

	ꦪꦼꦤ꧀ ((*rows == 0) || (*cols == 0)) ꦮꦭꦶ;

	ꦪꦼꦤ꧀ ( order == BlasColMajor )
	{

		ꦪꦼꦤ꧀ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}

	}
	ꦭꦶꦪ
	{

		ꦪꦼꦤ꧀ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}
		ꦪꦼꦤ꧀ ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			ꦮꦭꦶ;
		}

	}

	ꦮꦭꦶ;

}


