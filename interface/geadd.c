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

#もし defined(DOUBLE)
#define ERROR_NAME "DGEADD "
#否则
#define ERROR_NAME "SGEADD "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
                        FLOAT *BETA,  FLOAT *c, blasint *LDC, char *TRANS_A,char *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
char transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
char transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;
 もし(transa){
  もし (lda < MAX(1, n))info = 5;


  }	否则{

   もし (lda < MAX(1, m))	info = 5;
  }
 

  もし(transc){
        もし (ldc < MAX(1, n))	info = 8;


   }否则{
       
      もし (ldc < MAX(1, m))	info = 8;
          }

  もし (n < 0)		info = 2;
  もし (m < 0)		info = 1;

  もし (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    戻る;
  }

#否则
void CNAME( enum CBLAS_ORDER order,enum CBLAS_TRANSPOSE transa, enum CBLAS_TRANSPOSE transc,blasint m,  blasint n,  FLOAT alpha, FLOAT *a,  blasint lda, FLOAT beta, 
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

  もし (order == CblasColMajor) {

    info = -1;
もし ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      もし (ldc < MAX(1, m))  info = 8;
    } 否则 {
      もし (ldc < MAX(1, n))  info = 8;
    }

    もし ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      もし (lda < MAX(1, m))  info = 5;
    } 否则 {
      もし (lda < MAX(1, n))  info = 5;
    }
    もし (n < 0)		  info = 2;
    もし (m < 0)		  info = 1;

  }

  もし (order == CblasRowMajor) {
    info = -1;

    t = n;
    n = m;
    m = t;

もし ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      もし (ldc < MAX(1, m))  info = 8;
    } 否则 {
      もし (ldc < MAX(1, n))  info = 8;
    }

    もし ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      もし (lda < MAX(1, m))  info = 5;
    } 否则 {
      もし (lda < MAX(1, n))  info = 5;
    }
    もし (n < 0)		  info = 1;
    もし (m < 0)		  info = 2;
  }

  もし (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    戻る;
  }

#endif

  もし ((m==0) || (n==0)) 戻る;


  IDEBUG_START;

  FUNCTION_PROFILE_START();


  GEADD_K(m,n,alpha, a, lda, beta, c, ldc,
#ifdef CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#否则
          transa, transc
#endif
          ); 


  FUNCTION_PROFILE_END(1, 2* m * n ,  2 * m * n);

  IDEBUG_END;

  戻る;

}
