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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#যদি defined(DOUBLE)
#define ERROR_NAME "DGEADD "
#অন্যথা
#define ERROR_NAME "SGEADD "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
                        FLOAT *BETA,  FLOAT *c, blasint *LDC, অক্ষর *TRANS_A,অক্ষর *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
অক্ষর transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
অক্ষর transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;
 যদি(transa){
  যদি (lda < MAX(1, n))info = 5;


  }	অন্যথা{

   যদি (lda < MAX(1, m))	info = 5;
  }
 

  যদি(transc){
        যদি (ldc < MAX(1, n))	info = 8;


   }অন্যথা{
       
      যদি (ldc < MAX(1, m))	info = 8;
          }

  যদি (n < 0)		info = 2;
  যদি (m < 0)		info = 1;

  যদি (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#অন্যথা
শূন্য CNAME( গণনা CBLAS_ORDER order,গণনা CBLAS_TRANSPOSE transa, গণনা CBLAS_TRANSPOSE transc,blasint m,  blasint n,  FLOAT alpha, FLOAT *a,  blasint lda, FLOAT beta, 
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

  যদি (order == CblasColMajor) {

    info = -1;
যদি ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      যদি (ldc < MAX(1, m))  info = 8;
    } অন্যথা {
      যদি (ldc < MAX(1, n))  info = 8;
    }

    যদি ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      যদি (lda < MAX(1, m))  info = 5;
    } অন্যথা {
      যদি (lda < MAX(1, n))  info = 5;
    }
    যদি (n < 0)		  info = 2;
    যদি (m < 0)		  info = 1;

  }

  যদি (order == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;

যদি ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      যদি (ldc < MAX(1, m))  info = 8;
    } অন্যথা {
      যদি (ldc < MAX(1, n))  info = 8;
    }

    যদি ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      যদি (lda < MAX(1, m))  info = 5;
    } অন্যথা {
      যদি (lda < MAX(1, n))  info = 5;
    }
    যদি (n < 0)		  info = 1;
    যদি (m < 0)		  info = 2;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#endif

  যদি ((m==0) || (n==0)) ফেরত;


  IDEBUG_START;

  FUNCTION_PROFILE_START();


  GEADD_K(m,n,alpha, a, lda, beta, c, ldc,
#ifdef CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#অন্যথা
          transa, transc
#endif
          ); 


  FUNCTION_PROFILE_END(1, 2* m * n ,  2 * m * n);

  IDEBUG_END;

  ফেরত;

}
