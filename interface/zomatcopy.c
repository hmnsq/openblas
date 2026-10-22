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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#यदि defined(DOUBLE)
#define ERROR_NAME "ZOMATCOPY"
#अन्यथा
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
शून्य NAME( वर्ण* ORDER, वर्ण* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	वर्ण Order, Trans;
	पूर्णाङ्क order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	यदि ( Order == 'C' ) order = BlasColMajor;
	यदि ( Order == 'R' ) order = BlasRowMajor;
	यदि ( Trans == 'N' ) trans = BlasNoTrans;
	यदि ( Trans == 'T' ) trans = BlasTrans;
	यदि ( Trans == 'C' ) trans = BlasTransConj;
	यदि ( Trans == 'R' ) trans = BlasConj;

#अन्यथा 
शून्य CNAME(गणना CBLAS_ORDER CORDER, गणना CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	पूर्णाङ्क order=-1,trans=-1;
	blasint info = -1;

	यदि ( CORDER == CblasColMajor ) order = BlasColMajor; 
	यदि ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	यदि ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	यदि ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	यदि ( CTRANS == CblasTrans) trans = BlasTrans; 
	यदि ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif
	यदि ( order == BlasColMajor)
	{
        	यदि ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	यदि ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	यदि ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	यदि ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	यदि ( order == BlasRowMajor)
	{
        	यदि ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	यदि ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	यदि ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	यदि ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	यदि ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	यदि ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	यदि ( *cols < 0 ) info = 4;
	यदि ( *rows < 0 ) info = 3;
	यदि ( trans < 0 ) info = 2;
	यदि ( order < 0 ) info = 1;

	यदि (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    		प्रत्यागम;
  	}

	यदि ((*rows == 0) || (*cols == 0)) प्रत्यागम;

	यदि ( order == BlasColMajor )
	{

		यदि ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}

	}
	अन्यथा
	{

		यदि ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}
		यदि ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *ldb );
			प्रत्यागम;
		}

	}

	प्रत्यागम;

}


