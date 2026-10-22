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

#ضمن <stdio.h>
#ضمن "common.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#إذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGEADD "
#وإلا
#تعریف_کړه ERROR_NAME "SGEADD "
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
                        FLOAT *BETA,  FLOAT *c, blasint *LDC, حرف *TRANS_A,حرف *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
حرف transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
حرف transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;
 إذا(transa){
  إذا (lda < MAX(1, n))info = 5;


  }	وإلا{

   إذا (lda < MAX(1, m))	info = 5;
  }
 

  إذا(transc){
        إذا (ldc < MAX(1, n))	info = 8;


   }وإلا{
       
      إذا (ldc < MAX(1, m))	info = 8;
          }

  إذا (n < 0)		info = 2;
  إذا (m < 0)		info = 1;

  إذا (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا
فراغ CNAME( تعداد CBLAS_ترتیب ترتیب,تعداد CBLAS_TRANSPOSE transa, تعداد CBLAS_TRANSPOSE transc,blasint m,  blasint n,  FLOAT alpha, FLOAT *a,  blasint lda, FLOAT beta, 
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

  إذا (ترتیب == Cblas_ستون_Major) {

    info = -1;
إذا ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      إذا (ldc < MAX(1, m))  info = 8;
    } وإلا {
      إذا (ldc < MAX(1, n))  info = 8;
    }

    إذا ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      إذا (lda < MAX(1, m))  info = 5;
    } وإلا {
      إذا (lda < MAX(1, n))  info = 5;
    }
    إذا (n < 0)		  info = 2;
    إذا (m < 0)		  info = 1;

  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

إذا ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      إذا (ldc < MAX(1, m))  info = 8;
    } وإلا {
      إذا (ldc < MAX(1, n))  info = 8;
    }

    إذا ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      إذا (lda < MAX(1, m))  info = 5;
    } وإلا {
      إذا (lda < MAX(1, n))  info = 5;
    }
    إذا (n < 0)		  info = 1;
    إذا (m < 0)		  info = 2;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#د_که_پای

  إذا ((m==0) || (n==0)) أرجع;


  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();


  GEADD_K(m,n,alpha, a, lda, beta, c, ldc,
#که_تعریف_شوی CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#وإلا
          transa, transc
#د_که_پای
          ); 


  FUNCTION_PROFILE_پای_ټکی(1, 2* m * n ,  2 * m * n);

  IDEBUG_پای_ټکی;

  أرجع;

}
