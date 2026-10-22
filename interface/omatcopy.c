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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#എങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DOMATCOPY"
#അല്ലെങ്കിൽ
#define ERROR_NAME "SOMATCOPY"
#endif

#define BlasRowMajor 0
#define BlasColMajor 1
#define BlasNoTrans  0
#define BlasTrans    1

OPENBLAS_EXPORT

#ifndef CBLAS 
ശൂന്യം NAME( അക്ഷരം* ORDER, അക്ഷരം* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	അക്ഷരം Order, Trans;
	പൂർണ്ണസംഖ്യ order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	എങ്കിൽ ( Order == 'C' ) order = BlasColMajor;
	എങ്കിൽ ( Order == 'R' ) order = BlasRowMajor;
	എങ്കിൽ ( Trans == 'N' ) trans = BlasNoTrans;
	എങ്കിൽ ( Trans == 'R' ) trans = BlasNoTrans;
	എങ്കിൽ ( Trans == 'T' ) trans = BlasTrans;
	എങ്കിൽ ( Trans == 'C' ) trans = BlasTrans;
#അല്ലെങ്കിൽ 
ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER CORDER, എണ്ണൽ CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	FLOAT   *alpha; 
	പൂർണ്ണസംഖ്യ order=-1,trans=-1;
	blasint info = -1;

	എങ്കിൽ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	എങ്കിൽ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	എങ്കിൽ ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	എങ്കിൽ ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#endif
	എങ്കിൽ ( order == BlasColMajor)
	{
        	എങ്കിൽ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 9;
        	എങ്കിൽ ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	എങ്കിൽ ( order == BlasRowMajor)
	{
        	എങ്കിൽ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 9;
        	എങ്കിൽ ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	എങ്കിൽ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	എങ്കിൽ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	എങ്കിൽ ( *cols < 0 ) info = 4;
	എങ്കിൽ ( *rows < 0 ) info = 3;
	എങ്കിൽ ( trans < 0 ) info = 2;
	എങ്കിൽ ( order < 0 ) info = 1;

	എങ്കിൽ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    		തിരികെ;
  	}

	എങ്കിൽ ((*rows == 0) || (*cols == 0)) തിരികെ;

	എങ്കിൽ ( order == BlasColMajor )
	{
		എങ്കിൽ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		അല്ലെങ്കിൽ
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}
	അല്ലെങ്കിൽ
	{
		എങ്കിൽ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		അല്ലെങ്കിൽ
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}

	തിരികെ;

}


