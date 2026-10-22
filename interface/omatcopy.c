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
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DOMATCOPY"
#ورنہ
#تعریف_کرو ERROR_NAME "SOMATCOPY"
#اگر_ختم

#تعریف_کرو Blas_صف_Major 0
#تعریف_کرو Blas_ستون_Major 1
#تعریف_کرو BlasNoTrans  0
#تعریف_کرو BlasTrans    1

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS 
خالی NAME( حرف* ترتیب_2, حرف* TRANS, blasint *صفیں, blasint *ستون_تعداد, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	حرف ترتیب_3, Trans;
	صحیح ترتیب=-1,trans=-1;
	blasint info = -1;

	ترتیب_3 = *ترتیب_2;
	Trans = *TRANS;

	TOUPPER(ترتیب_3);
	TOUPPER(Trans);

	اگر ( ترتیب_3 == 'C' ) ترتیب = Blas_ستون_Major;
	اگر ( ترتیب_3 == 'R' ) ترتیب = Blas_صف_Major;
	اگر ( Trans == 'N' ) trans = BlasNoTrans;
	اگر ( Trans == 'R' ) trans = BlasNoTrans;
	اگر ( Trans == 'T' ) trans = BlasTrans;
	اگر ( Trans == 'C' ) trans = BlasTrans;
#ورنہ 
خالی CNAME(شمار CBLAS_ترتیب CORDER, شمار CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *صفیں, *ستون_تعداد, *lda, *ldb; 
	FLOAT   *alpha; 
	صحیح ترتیب=-1,trans=-1;
	blasint info = -1;

	اگر ( CORDER == Cblas_ستون_Major ) ترتیب = Blas_ستون_Major; 
	اگر ( CORDER == Cblas_صف_Major ) ترتیب = Blas_صف_Major; 

	اگر ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	اگر ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	صفیں = &crows; 
	ستون_تعداد = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#اگر_ختم
	اگر ( ترتیب == Blas_ستون_Major)
	{
        	اگر ( trans == BlasNoTrans  &&  *ldb < MAX(1,*صفیں) ) info = 9;
        	اگر ( trans == BlasTrans    &&  *ldb < MAX(1,*ستون_تعداد) ) info = 9;
	}
	اگر ( ترتیب == Blas_صف_Major)
	{
        	اگر ( trans == BlasNoTrans  &&  *ldb < MAX(1,*ستون_تعداد) ) info = 9;
        	اگر ( trans == BlasTrans    &&  *ldb < MAX(1,*صفیں) ) info = 9;
	}

	اگر ( ترتیب == Blas_ستون_Major &&  *lda < MAX(1,*صفیں) ) info = 7;
	اگر ( ترتیب == Blas_صف_Major &&  *lda < MAX(1,*ستون_تعداد) ) info = 7;
	اگر ( *ستون_تعداد < 0 ) info = 4;
	اگر ( *صفیں < 0 ) info = 3;
	اگر ( trans < 0 ) info = 2;
	اگر ( ترتیب < 0 ) info = 1;

	اگر (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    		واپس;
  	}

	اگر ((*صفیں == 0) || (*ستون_تعداد == 0)) واپس;

	اگر ( ترتیب == Blas_ستون_Major )
	{
		اگر ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*صفیں, *ستون_تعداد, *alpha, a, *lda, b, *ldb );
		}
		ورنہ
		{
			OMATCOPY_K_CT(*صفیں, *ستون_تعداد, *alpha, a, *lda, b, *ldb );
		}
	}
	ورنہ
	{
		اگر ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*صفیں, *ستون_تعداد, *alpha, a, *lda, b, *ldb );
		}
		ورنہ
		{
			OMATCOPY_K_RT(*صفیں, *ستون_تعداد, *alpha, a, *lda, b, *ldb );
		}
	}

	واپس;

}


