/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEADD "
#अन्यथा
#परिभाषय ERROR_NAME "SGEADD "
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
                        FLOAT *BETA,  FLOAT *c, blasint *LDC, वर्ण *TRANS_A,वर्ण *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
वर्ण transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
वर्ण transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;
 यदि(transa){
  यदि (lda < MAX(1, n))info = 5;


  }	अन्यथा{

   यदि (lda < MAX(1, m))	info = 5;
  }
 

  यदि(transc){
        यदि (ldc < MAX(1, n))	info = 8;


   }अन्यथा{
       
      यदि (ldc < MAX(1, m))	info = 8;
          }

  यदि (n < 0)		info = 2;
  यदि (m < 0)		info = 1;

  यदि (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा
शून्य CNAME( गणना CBLAS_क्रमः क्रमः,गणना CBLAS_TRANSPOSE transa, गणना CBLAS_TRANSPOSE transc,blasint m,  blasint n,  FLOAT alpha, FLOAT *a,  blasint lda, FLOAT beta, 
		  FLOAT *c,  blasint ldc)
{
/* 
void CNAME(enum CBLAS_ORDER order,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT beta,
	   FLOAT  *c, blasint ldc){ */

  blasint info, t;

  PRINT_DEBUG_CNAME;

  info  =  0;

  यदि (क्रमः == CblasColMajor) {

    info = -1;
यदि ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      यदि (ldc < MAX(1, m))  info = 8;
    } अन्यथा {
      यदि (ldc < MAX(1, n))  info = 8;
    }

    यदि ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      यदि (lda < MAX(1, m))  info = 5;
    } अन्यथा {
      यदि (lda < MAX(1, n))  info = 5;
    }
    यदि (n < 0)		  info = 2;
    यदि (m < 0)		  info = 1;

  }

  यदि (क्रमः == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;

यदि ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      यदि (ldc < MAX(1, m))  info = 8;
    } अन्यथा {
      यदि (ldc < MAX(1, n))  info = 8;
    }

    यदि ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      यदि (lda < MAX(1, m))  info = 5;
    } अन्यथा {
      यदि (lda < MAX(1, n))  info = 5;
    }
    यदि (n < 0)		  info = 1;
    यदि (m < 0)		  info = 2;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  यदि ((m==0) || (n==0)) प्रत्यागम;


  IDEBUG_START;

  FUNCTION_PROFILE_START();


  GEADD_K(m,n,alpha, a, lda, beta, c, ldc,
#यदि_परिभाषितम् CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#अन्यथा
          transa, transc
#यदि_अन्तः
          ); 


  FUNCTION_PROFILE_अन्तः(1, 2* m * n ,  2 * m * n);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
