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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#إذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DOMATCOPY"
#وإلا
#تعریف_کړه ERROR_NAME "SOMATCOPY"
#د_که_پای

#تعریف_کړه Blas_لیکه_Major 0
#تعریف_کړه Blas_ستون_Major 1
#تعریف_کړه BlasNoTrans  0
#تعریف_کړه BlasTrans    1

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS 
فراغ NAME( حرف* ترتیب_2, حرف* TRANS, blasint *لیکې, blasint *ستون_شمېر, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	حرف ترتیب_3, Trans;
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	ترتیب_3 = *ترتیب_2;
	Trans = *TRANS;

	TOUPPER(ترتیب_3);
	TOUPPER(Trans);

	إذا ( ترتیب_3 == 'C' ) ترتیب = Blas_ستون_Major;
	إذا ( ترتیب_3 == 'R' ) ترتیب = Blas_لیکه_Major;
	إذا ( Trans == 'N' ) trans = BlasNoTrans;
	إذا ( Trans == 'R' ) trans = BlasNoTrans;
	إذا ( Trans == 'T' ) trans = BlasTrans;
	إذا ( Trans == 'C' ) trans = BlasTrans;
#وإلا 
فراغ CNAME(تعداد CBLAS_ترتیب CORDER, تعداد CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *لیکې, *ستون_شمېر, *lda, *ldb; 
	FLOAT   *alpha; 
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	إذا ( CORDER == Cblas_ستون_Major ) ترتیب = Blas_ستون_Major; 
	إذا ( CORDER == Cblas_لیکه_Major ) ترتیب = Blas_لیکه_Major; 

	إذا ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	إذا ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	لیکې = &crows; 
	ستون_شمېر = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#د_که_پای
	إذا ( ترتیب == Blas_ستون_Major)
	{
        	إذا ( trans == BlasNoTrans  &&  *ldb < MAX(1,*لیکې) ) info = 9;
        	إذا ( trans == BlasTrans    &&  *ldb < MAX(1,*ستون_شمېر) ) info = 9;
	}
	إذا ( ترتیب == Blas_لیکه_Major)
	{
        	إذا ( trans == BlasNoTrans  &&  *ldb < MAX(1,*ستون_شمېر) ) info = 9;
        	إذا ( trans == BlasTrans    &&  *ldb < MAX(1,*لیکې) ) info = 9;
	}

	إذا ( ترتیب == Blas_ستون_Major &&  *lda < MAX(1,*لیکې) ) info = 7;
	إذا ( ترتیب == Blas_لیکه_Major &&  *lda < MAX(1,*ستون_شمېر) ) info = 7;
	إذا ( *ستون_شمېر < 0 ) info = 4;
	إذا ( *لیکې < 0 ) info = 3;
	إذا ( trans < 0 ) info = 2;
	إذا ( ترتیب < 0 ) info = 1;

	إذا (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    		أرجع;
  	}

	إذا ((*لیکې == 0) || (*ستون_شمېر == 0)) أرجع;

	إذا ( ترتیب == Blas_ستون_Major )
	{
		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*لیکې, *ستون_شمېر, *alpha, a, *lda, b, *ldb );
		}
		وإلا
		{
			OMATCOPY_K_CT(*لیکې, *ستون_شمېر, *alpha, a, *lda, b, *ldb );
		}
	}
	وإلا
	{
		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*لیکې, *ستون_شمېر, *alpha, a, *lda, b, *ldb );
		}
		وإلا
		{
			OMATCOPY_K_RT(*لیکې, *ستون_شمېر, *alpha, a, *lda, b, *ldb );
		}
	}

	أرجع;

}


