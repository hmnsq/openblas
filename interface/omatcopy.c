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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#إذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DOMATCOPY"
#وإلا
#تعریف_کن ERROR_NAME "SOMATCOPY"
#پایان_اگر

#تعریف_کن Blas_ردیف_Major 0
#تعریف_کن Blas_ستون_Major 1
#تعریف_کن BlasNoTrans  0
#تعریف_کن BlasTrans    1

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS 
فراغ NAME( حرف* ترتیب_2, حرف* TRANS, blasint *ردیفها, blasint *تعداد_ستون, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	حرف ترتیب_3, Trans;
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	ترتیب_3 = *ترتیب_2;
	Trans = *TRANS;

	TOUPPER(ترتیب_3);
	TOUPPER(Trans);

	إذا ( ترتیب_3 == 'C' ) ترتیب = Blas_ستون_Major;
	إذا ( ترتیب_3 == 'R' ) ترتیب = Blas_ردیف_Major;
	إذا ( Trans == 'N' ) trans = BlasNoTrans;
	إذا ( Trans == 'R' ) trans = BlasNoTrans;
	إذا ( Trans == 'T' ) trans = BlasTrans;
	إذا ( Trans == 'C' ) trans = BlasTrans;
#وإلا 
فراغ CNAME(تعداد CBLAS_ترتیب CORDER, تعداد CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *ردیفها, *تعداد_ستون, *lda, *ldb; 
	FLOAT   *alpha; 
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	إذا ( CORDER == Cblas_ستون_Major ) ترتیب = Blas_ستون_Major; 
	إذا ( CORDER == Cblas_ردیف_Major ) ترتیب = Blas_ردیف_Major; 

	إذا ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	إذا ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	ردیفها = &crows; 
	تعداد_ستون = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#پایان_اگر
	إذا ( ترتیب == Blas_ستون_Major)
	{
        	إذا ( trans == BlasNoTrans  &&  *ldb < MAX(1,*ردیفها) ) info = 9;
        	إذا ( trans == BlasTrans    &&  *ldb < MAX(1,*تعداد_ستون) ) info = 9;
	}
	إذا ( ترتیب == Blas_ردیف_Major)
	{
        	إذا ( trans == BlasNoTrans  &&  *ldb < MAX(1,*تعداد_ستون) ) info = 9;
        	إذا ( trans == BlasTrans    &&  *ldb < MAX(1,*ردیفها) ) info = 9;
	}

	إذا ( ترتیب == Blas_ستون_Major &&  *lda < MAX(1,*ردیفها) ) info = 7;
	إذا ( ترتیب == Blas_ردیف_Major &&  *lda < MAX(1,*تعداد_ستون) ) info = 7;
	إذا ( *تعداد_ستون < 0 ) info = 4;
	إذا ( *ردیفها < 0 ) info = 3;
	إذا ( trans < 0 ) info = 2;
	إذا ( ترتیب < 0 ) info = 1;

	إذا (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    		أرجع;
  	}

	إذا ((*ردیفها == 0) || (*تعداد_ستون == 0)) أرجع;

	إذا ( ترتیب == Blas_ستون_Major )
	{
		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*ردیفها, *تعداد_ستون, *alpha, a, *lda, b, *ldb );
		}
		وإلا
		{
			OMATCOPY_K_CT(*ردیفها, *تعداد_ستون, *alpha, a, *lda, b, *ldb );
		}
	}
	وإلا
	{
		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*ردیفها, *تعداد_ستون, *alpha, a, *lda, b, *ldb );
		}
		وإلا
		{
			OMATCOPY_K_RT(*ردیفها, *تعداد_ستون, *alpha, a, *lda, b, *ldb );
		}
	}

	أرجع;

}


