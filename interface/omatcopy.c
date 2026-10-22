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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdlib.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DOMATCOPY"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SOMATCOPY"
#ਜੇ_ਅੰਤ

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ Blas_ਪੰਕਤੀ_Major 0
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ Blas_ਖਾਨਾ_Major 1
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ BlasNoTrans  0
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ BlasTrans    1

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS 
ਖਾਲੀ NAME( ਅੱਖਰ* ਕ੍ਰਮ_2, ਅੱਖਰ* TRANS, blasint *ਪੰਕਤੀਆਂ, blasint *ਖਾਨੇ_ਗਿਣਤੀ, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	ਅੱਖਰ ਕ੍ਰਮ_3, Trans;
	ਪੂਰਨ ਕ੍ਰਮ=-1,trans=-1;
	blasint info = -1;

	ਕ੍ਰਮ_3 = *ਕ੍ਰਮ_2;
	Trans = *TRANS;

	TOUPPER(ਕ੍ਰਮ_3);
	TOUPPER(Trans);

	ਜੇ ( ਕ੍ਰਮ_3 == 'C' ) ਕ੍ਰਮ = Blas_ਖਾਨਾ_Major;
	ਜੇ ( ਕ੍ਰਮ_3 == 'R' ) ਕ੍ਰਮ = Blas_ਪੰਕਤੀ_Major;
	ਜੇ ( Trans == 'N' ) trans = BlasNoTrans;
	ਜੇ ( Trans == 'R' ) trans = BlasNoTrans;
	ਜੇ ( Trans == 'T' ) trans = BlasTrans;
	ਜੇ ( Trans == 'C' ) trans = BlasTrans;
#ਨਹੀਂਤਾਂ 
ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ CORDER, ਗਿਣਤੀ CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *ਪੰਕਤੀਆਂ, *ਖਾਨੇ_ਗਿਣਤੀ, *lda, *ldb; 
	FLOAT   *alpha; 
	ਪੂਰਨ ਕ੍ਰਮ=-1,trans=-1;
	blasint info = -1;

	ਜੇ ( CORDER == Cblas_ਖਾਨਾ_Major ) ਕ੍ਰਮ = Blas_ਖਾਨਾ_Major; 
	ਜੇ ( CORDER == Cblas_ਪੰਕਤੀ_Major ) ਕ੍ਰਮ = Blas_ਪੰਕਤੀ_Major; 

	ਜੇ ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	ਜੇ ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	ਪੰਕਤੀਆਂ = &crows; 
	ਖਾਨੇ_ਗਿਣਤੀ = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#ਜੇ_ਅੰਤ
	ਜੇ ( ਕ੍ਰਮ == Blas_ਖਾਨਾ_Major)
	{
        	ਜੇ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*ਪੰਕਤੀਆਂ) ) info = 9;
        	ਜੇ ( trans == BlasTrans    &&  *ldb < MAX(1,*ਖਾਨੇ_ਗਿਣਤੀ) ) info = 9;
	}
	ਜੇ ( ਕ੍ਰਮ == Blas_ਪੰਕਤੀ_Major)
	{
        	ਜੇ ( trans == BlasNoTrans  &&  *ldb < MAX(1,*ਖਾਨੇ_ਗਿਣਤੀ) ) info = 9;
        	ਜੇ ( trans == BlasTrans    &&  *ldb < MAX(1,*ਪੰਕਤੀਆਂ) ) info = 9;
	}

	ਜੇ ( ਕ੍ਰਮ == Blas_ਖਾਨਾ_Major &&  *lda < MAX(1,*ਪੰਕਤੀਆਂ) ) info = 7;
	ਜੇ ( ਕ੍ਰਮ == Blas_ਪੰਕਤੀ_Major &&  *lda < MAX(1,*ਖਾਨੇ_ਗਿਣਤੀ) ) info = 7;
	ਜੇ ( *ਖਾਨੇ_ਗਿਣਤੀ < 0 ) info = 4;
	ਜੇ ( *ਪੰਕਤੀਆਂ < 0 ) info = 3;
	ਜੇ ( trans < 0 ) info = 2;
	ਜੇ ( ਕ੍ਰਮ < 0 ) info = 1;

	ਜੇ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    		ਵਾਪਸ;
  	}

	ਜੇ ((*ਪੰਕਤੀਆਂ == 0) || (*ਖਾਨੇ_ਗਿਣਤੀ == 0)) ਵਾਪਸ;

	ਜੇ ( ਕ੍ਰਮ == Blas_ਖਾਨਾ_Major )
	{
		ਜੇ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*ਪੰਕਤੀਆਂ, *ਖਾਨੇ_ਗਿਣਤੀ, *alpha, a, *lda, b, *ldb );
		}
		ਨਹੀਂਤਾਂ
		{
			OMATCOPY_K_CT(*ਪੰਕਤੀਆਂ, *ਖਾਨੇ_ਗਿਣਤੀ, *alpha, a, *lda, b, *ldb );
		}
	}
	ਨਹੀਂਤਾਂ
	{
		ਜੇ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*ਪੰਕਤੀਆਂ, *ਖਾਨੇ_ਗਿਣਤੀ, *alpha, a, *lda, b, *ldb );
		}
		ਨਹੀਂਤਾਂ
		{
			OMATCOPY_K_RT(*ਪੰਕਤੀਆਂ, *ਖਾਨੇ_ਗਿਣਤੀ, *alpha, a, *lda, b, *ldb );
		}
	}

	ਵਾਪਸ;

}


