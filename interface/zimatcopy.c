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
 * 2014-06-10 Saar
 * 2015-09-07 grisuthedragon 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ຖ້າ defined(DOUBLE)
#define ERROR_NAME "ZIMATCOPY"
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "CIMATCOPY"
#endif

#define BlasRowMajor     0
#define BlasColMajor     1
#define BlasNoTrans      0
#define BlasTrans        1
#define BlasTransConj    2
#define BlasConj         3

#define NEW_IMATCOPY 

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME( var_na* ORDER, var_na* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, blasint *ldb)
{

	var_na Order, Trans;
	ຈຳນວນເຕັມ order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	ຖ້າ ( Order == 'C' ) order = BlasColMajor;
	ຖ້າ ( Order == 'R' ) order = BlasRowMajor;
	ຖ້າ ( Trans == 'N' ) trans = BlasNoTrans;
	ຖ້າ ( Trans == 'T' ) trans = BlasTrans;
	ຖ້າ ( Trans == 'C' ) trans = BlasTransConj;
	ຖ້າ ( Trans == 'R' ) trans = BlasConj;

#ບໍ່ດັ່ງນັ້ນ 
void CNAME( ga_nanaa CBLAS_ORDER CORDER, ga_nanaa CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT *alpha, FLOAT *a, blasint clda, blasint cldb)
{

	blasint *rows, *cols, *lda, *ldb; 
	ຈຳນວນເຕັມ order=-1,trans=-1;
	blasint info = -1;
	FLOAT *b;
	size_t msize;

	ຖ້າ ( CORDER == CblasColMajor ) order = BlasColMajor; 
	ຖ້າ ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	ຖ້າ ( CTRANS == CblasNoTrans) trans = BlasNoTrans; 
	ຖ້າ ( CTRANS == CblasConjNoTrans ) trans = BlasConj; 
	ຖ້າ ( CTRANS == CblasTrans) trans = BlasTrans; 
	ຖ້າ ( CTRANS == CblasConjTrans) trans = BlasTransConj; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
#endif

	ຖ້າ ( order == BlasColMajor)
	{
        	ຖ້າ ( trans == BlasNoTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ຖ້າ ( trans == BlasConj         &&  *ldb < MAX(1,*rows) ) info = 9;
        	ຖ້າ ( trans == BlasTrans        &&  *ldb < MAX(1,*cols) ) info = 9;
        	ຖ້າ ( trans == BlasTransConj    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	ຖ້າ ( order == BlasRowMajor)
	{
        	ຖ້າ ( trans == BlasNoTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
        	ຖ້າ ( trans == BlasConj       &&  *ldb < MAX(1,*cols) ) info = 9;
        	ຖ້າ ( trans == BlasTrans      &&  *ldb < MAX(1,*rows) ) info = 9;
        	ຖ້າ ( trans == BlasTransConj  &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	ຖ້າ ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	ຖ້າ ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	ຖ້າ ( *cols < 0 ) info = 4;
	ຖ້າ ( *rows < 0 ) info = 3;
	ຖ້າ ( trans < 0 ) info = 2;
	ຖ້າ ( order < 0 ) info = 1;

	ຖ້າ (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    		ສົ່ງຄືນ;
  	}

	ຖ້າ ((*rows == 0) || (*cols == 0)) ສົ່ງຄືນ;

#ifdef NEW_IMATCOPY
    ຖ້າ (*lda == *ldb ) {
        ຖ້າ ( order == BlasColMajor )
        {

            ຖ້າ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasConj )
            {
                IMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasTrans && *rows == *cols )
            {
                IMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasTransConj && *rows == *cols )
            {
                IMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }

        }
        ບໍ່ດັ່ງນັ້ນ
        {

            ຖ້າ ( trans == BlasNoTrans )
            {
                IMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasConj )
            {
                IMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasTrans && *rows == *cols )
            {
                IMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }
            ຖ້າ ( trans == BlasTransConj && *rows == *cols )
            {
                IMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda );
                ສົ່ງຄືນ;
            }

        }
    }
#endif

		ຖ້າ ( *rows >  *cols )
                msize = (size_t)(*rows) * (*ldb)  * parimaa_na(FLOAT) * 2;
        ບໍ່ດັ່ງນັ້ນ
                msize = (size_t)(*cols) * (*ldb)  * parimaa_na(FLOAT) * 2;

	b = malloc(msize);
	ຖ້າ ( b == NULL )
	{
		printf("Memory alloc failed in zimatcopy\n");
		exit(1);
	}

	ຖ້າ ( order == BlasColMajor )
	{

		ຖ້າ ( trans == BlasNoTrans )
		{
	  		OMATCOPY_K_CN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasConj )
		{
			OMATCOPY_K_CNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_CN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasTrans )
		{
			OMATCOPY_K_CT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasTransConj )
		{
			OMATCOPY_K_CTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_CN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}

	}
	ບໍ່ດັ່ງນັ້ນ
	{

		ຖ້າ ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasConj )
		{
			OMATCOPY_K_RNC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *cols );
	  		OMATCOPY_K_RN(*rows, *cols, (FLOAT) 1.0, (FLOAT) 0.0 , b, *cols, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasTrans )
		{
			OMATCOPY_K_RT(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}
		ບໍ່ດັ່ງນັ້ນ ຖ້າ ( trans == BlasTransConj )
		{
			OMATCOPY_K_RTC(*rows, *cols, alpha[0], alpha[1], a, *lda, b, *rows );
	  		OMATCOPY_K_RN(*cols, *rows, (FLOAT) 1.0, (FLOAT) 0.0 , b, *rows, a, *ldb );
		}

	}

	free(b);
	ສົ່ງຄືນ;

}


