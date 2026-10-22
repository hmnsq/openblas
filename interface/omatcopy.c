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

#အကယ်၍ defined(DOUBLE)
#define ERROR_NAME "DOMATCOPY"
#မဟုတ်ပါက
#define ERROR_NAME "SOMATCOPY"
#endif

#define BlasRowMajor 0
#define BlasColMajor 1
#define BlasNoTrans  0
#define BlasTrans    1

OPENBLAS_EXPORT

#ifndef CBLAS 
void NAME( var_na* ORDER, var_na* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	var_na Order, Trans;
	ကိန်းပြည့် order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	အကယ်၍ ( Order == 'C' ) order = BlasColMajor;
	အကယ်၍ ( Order == 'R' ) order = BlasRowMajor;
	အကယ်၍ ( Trans == 'N' ) trans = BlasNoTrans;
	အကယ်၍ ( Trans == 'R' ) trans = BlasNoTrans;
	အကယ်၍ ( Trans == 'T' ) trans = BlasTrans;
	အကယ်၍ ( Trans == 'C' ) trans = BlasTrans;
#မဟုတ်ပါက 
void CNAME(ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	FLOAT   *alpha; 
	ကိန်းပြည့် order=-1,trans=-1;
	blasint info = -1;

	အကယ်၍ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	အကယ်၍ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	အကယ်၍ ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	အကယ်၍ ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#endif
	အကယ်၍ ( order == BlasColMajor)
	{
        	အကယ်၍ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 9;
        	အကယ်၍ ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	အကယ်၍ ( order == BlasRowMajor)
	{
        	အကယ်၍ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 9;
        	အကယ်၍ ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	အကယ်၍ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	အကယ်၍ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	အကယ်၍ ( *cols < 0 ) info = 4;
	အကယ်၍ ( *rows < 0 ) info = 3;
	အကယ်၍ ( trans < 0 ) info = 2;
	အကယ်၍ ( order < 0 ) info = 1;

	အကယ်၍ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ပြန်ပေး;
  	}

	အကယ်၍ ((*rows == 0) || (*cols == 0)) ပြန်ပေး;

	အကယ်၍ ( order == BlasColMajor )
	{
		အကယ်၍ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		မဟုတ်ပါက
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}
	မဟုတ်ပါက
	{
		အကယ်၍ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		မဟုတ်ပါက
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}

	ပြန်ပေး;

}


