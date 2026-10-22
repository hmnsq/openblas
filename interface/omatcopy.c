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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DOMATCOPY"
#अन्यथा
#घोषणा ERROR_NAME "SOMATCOPY"
#पूर्ण_अगर

#घोषणा BlasRowMajor 0
#घोषणा BlasColMajor 1
#घोषणा BlasNoTrans  0
#घोषणा BlasTrans    1

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS 
व्योम NAME( अक्षर* ORDER, अक्षर* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	अक्षर Order, Trans;
	पूर्णांक order=-1,trans=-1;
	blasint info = -1;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	अगर ( Order == 'C' ) order = BlasColMajor;
	अगर ( Order == 'R' ) order = BlasRowMajor;
	अगर ( Trans == 'N' ) trans = BlasNoTrans;
	अगर ( Trans == 'R' ) trans = BlasNoTrans;
	अगर ( Trans == 'T' ) trans = BlasTrans;
	अगर ( Trans == 'C' ) trans = BlasTrans;
#अन्यथा 
व्योम CNAME(क्रमागत CBLAS_ORDER CORDER, क्रमागत CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	FLOAT   *alpha; 
	पूर्णांक order=-1,trans=-1;
	blasint info = -1;

	अगर ( CORDER == CblasColMajor ) order = BlasColMajor; 
	अगर ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	अगर ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	अगर ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#पूर्ण_अगर
	अगर ( order == BlasColMajor)
	{
        	अगर ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 9;
        	अगर ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	अगर ( order == BlasRowMajor)
	{
        	अगर ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 9;
        	अगर ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	अगर ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	अगर ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	अगर ( *cols < 0 ) info = 4;
	अगर ( *rows < 0 ) info = 3;
	अगर ( trans < 0 ) info = 2;
	अगर ( order < 0 ) info = 1;

	अगर (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    		वापस;
  	}

	अगर ((*rows == 0) || (*cols == 0)) वापस;

	अगर ( order == BlasColMajor )
	{
		अगर ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		अन्यथा
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}
	अन्यथा
	{
		अगर ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		अन्यथा
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}

	वापस;

}


