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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ပေႃး defined(DOUBLE)
#define ERROR_NAME "ZGEADD "
#ပေႃးဢမ်ႇ
#define ERROR_NAME "CGEADD "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
		                  FLOAT *BETA,  FLOAT *c, blasint *LDC,var_na *TRANS_A,var_na *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 

  blasint info;
var_na transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
var_na transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;


ပေႃး(transa){
  ပေႃး (lda < MAX(1, n))info = 5;


  }	ပေႃးဢမ်ႇ{

   ပေႃး (lda < MAX(1, m))	info = 5;
  }
 

  ပေႃး(transc){
        ပေႃး (ldc < MAX(1, n))	info = 8;


   }ပေႃးဢမ်ႇ{
       
      ပေႃး (ldc < MAX(1, m))	info = 8;
          }


  ပေႃး (n < 0)		info = 2;
  ပေႃး (m < 0)		info = 1;

  ပေႃး (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

#ပေႃးဢမ်ႇ
void CNAME( ga_nanaa CBLAS_ORDER order,ga_nanaa CBLAS_TRANSPOSE transa, ga_nanaa CBLAS_TRANSPOSE transc, blasint m,  blasint n,  FLOAT *ALPHA, FLOAT *a,  blasint lda, FLOAT *BETA, 
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

  ပေႃး (order == CblasColMajor) {

    info = -1;

ပေႃး ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ပေႃး (ldc < MAX(1, m))  info = 8;
    } ပေႃးဢမ်ႇ {
      ပေႃး (ldc < MAX(1, n))  info = 8;
    }

    ပေႃး ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ပေႃး (lda < MAX(1, m))  info = 5;
    } ပေႃးဢမ်ႇ {
      ပေႃး (lda < MAX(1, n))  info = 5;
    }
    ပေႃး (n < 0)		  info = 2;
    ပေႃး (m < 0)		  info = 1;  
  }

  ပေႃး (order == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;
ပေႃး ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ပေႃး (ldc < MAX(1, m))  info = 8;
    } ပေႃးဢမ်ႇ {
      ပေႃး (ldc < MAX(1, n))  info = 8;
    }

    ပေႃး ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ပေႃး (lda < MAX(1, m))  info = 5;
    } ပေႃးဢမ်ႇ {
      ပေႃး (lda < MAX(1, n))  info = 5;
    }
    ပေႃး (n < 0)		  info = 1;
    ပေႃး (m < 0)		  info = 2;
  }

 
  

  ပေႃး (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

#endif

  ပေႃး ((m==0) || (n==0)) ပွၵ်ႈ;


  IDEBUG_START;

  FUNCTION_PROFILE_START();


  GEADD_K(m,n,ALPHA[0],ALPHA[1], a, lda, BETA[0], BETA[1], c, ldc,
        
          #ifdef CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#ပေႃးဢမ်ႇ
          transa, transc
#endif
          );


           


  FUNCTION_PROFILE_END(1, 2* m * n ,  2 * m * n);

  IDEBUG_END;

  ပွၵ်ႈ;

}
