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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#إذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZOMATCOPY"
#وإلا
#تعریف_کرو ERROR_NAME "COMATCOPY"
#جے_انت

#تعریف_کرو Blas_صف_Major     0
#تعریف_کرو Blas_خانہ_Major     1
#تعریف_کرو BlasNoTrans      0
#تعریف_کرو BlasTrans        1
#تعریف_کرو BlasTransConj    2
#تعریف_کرو BlasConj         3

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS
فراغ NAME( حرف* ترتیب_2, حرف* TRANS, blasint *صفاں, blasint *خانے_گنتی, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	حرف ترتیب_3, Trans;
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	ترتیب_3 = *ترتیب_2;
	Trans = *TRANS;

	TOUPPER(ترتیب_3);
	TOUPPER(Trans);

	إذا ( ترتیب_3 == 'C' ) ترتیب = Blas_خانہ_Major;
	إذا ( ترتیب_3 == 'R' ) ترتیب = Blas_صف_Major;
	إذا ( Trans == 'N' ) trans = BlasNoTrans;
	إذا ( Trans == 'T' ) trans = BlasTrans;
	إذا ( Trans == 'C' ) trans = BlasTransConj;
	إذا ( Trans == 'R' ) trans = BlasConj;

#وإلا 
فراغ CNAME(تعداد CBLAS_ترتیب CORDER, تعداد CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT  *alpha, FLOAT *a, blasint clda, FLOAT*b, blasint cldb)
{
	blasint *صفاں, *خانے_گنتی, *lda, *ldb; 
	صحيح ترتیب=-1,trans=-1;
	blasint info = -1;

	إذا ( CORDER == Cblas_خانہ_Major ) ترتیب = Blas_خانہ_Major; 
	إذا ( CORDER == Cblas_صف_Major ) ترتیب = Blas_صف_Major; 

	إذا ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	إذا ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	إذا ( CTRANS == CblasTrans) trans = BlasTrans; 
	إذا ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	صفاں = &crows; 
	خانے_گنتی = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#جے_انت
	إذا ( ترتیب == Blas_خانہ_Major)
	{
        	إذا ( trans == BlasNoTrans      &&  *ldb < MAX(1,*صفاں) ) info = 9;
        	إذا ( trans == BlasConj         &&  *ldb < MAX(1,*صفاں) ) info = 9;
        	إذا ( trans == BlasTrans        &&  *ldb < MAX(1,*خانے_گنتی) ) info = 9;
        	إذا ( trans == BlasTransConj    &&  *ldb < MAX(1,*خانے_گنتی) ) info = 9;
	}
	إذا ( ترتیب == Blas_صف_Major)
	{
        	إذا ( trans == BlasNoTrans    &&  *ldb < MAX(1,*خانے_گنتی) ) info = 9;
        	إذا ( trans == BlasConj       &&  *ldb < MAX(1,*خانے_گنتی) ) info = 9;
        	إذا ( trans == BlasTrans      &&  *ldb < MAX(1,*صفاں) ) info = 9;
        	إذا ( trans == BlasTransConj  &&  *ldb < MAX(1,*صفاں) ) info = 9;
	}

	إذا ( ترتیب == Blas_خانہ_Major &&  *lda < MAX(1,*صفاں) ) info = 7;
	إذا ( ترتیب == Blas_صف_Major &&  *lda < MAX(1,*خانے_گنتی) ) info = 7;
	إذا ( *خانے_گنتی < 0 ) info = 4;
	إذا ( *صفاں < 0 ) info = 3;
	إذا ( trans < 0 ) info = 2;
	إذا ( ترتیب < 0 ) info = 1;

	إذا (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    		أرجع;
  	}

	إذا ((*صفاں == 0) || (*خانے_گنتی == 0)) أرجع;

	إذا ( ترتیب == Blas_خانہ_Major )
	{

		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}

	}
	وإلا
	{

		إذا ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}
		إذا ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*صفاں, *خانے_گنتی, alpha[0], alpha[1], a, *lda, b, *ldb );
			أرجع;
		}

	}

	أرجع;

}


