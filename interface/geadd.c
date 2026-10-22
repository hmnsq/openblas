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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#ఒకవేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEADD "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEADD "
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
                        FLOAT *BETA,  FLOAT *c, blasint *LDC, అక్షరం *TRANS_A,అక్షరం *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
అక్షరం transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
అక్షరం transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;
 ఒకవేళ(transa){
  ఒకవేళ (lda < MAX(1, n))info = 5;


  }	లేకపోతే{

   ఒకవేళ (lda < MAX(1, m))	info = 5;
  }
 

  ఒకవేళ(transc){
        ఒకవేళ (ldc < MAX(1, n))	info = 8;


   }లేకపోతే{
       
      ఒకవేళ (ldc < MAX(1, m))	info = 8;
          }

  ఒకవేళ (n < 0)		info = 2;
  ఒకవేళ (m < 0)		info = 1;

  ఒకవేళ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే
శూన్య CNAME( గణన CBLAS_క్రమం క్రమం,గణన CBLAS_TRANSPOSE transa, గణన CBLAS_TRANSPOSE transc,blasint m,  blasint n,  FLOAT alpha, FLOAT *a,  blasint lda, FLOAT beta, 
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

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {

    info = -1;
ఒకవేళ ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ఒకవేళ (ldc < MAX(1, m))  info = 8;
    } లేకపోతే {
      ఒకవేళ (ldc < MAX(1, n))  info = 8;
    }

    ఒకవేళ ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ఒకవేళ (lda < MAX(1, m))  info = 5;
    } లేకపోతే {
      ఒకవేళ (lda < MAX(1, n))  info = 5;
    }
    ఒకవేళ (n < 0)		  info = 2;
    ఒకవేళ (m < 0)		  info = 1;

  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;

ఒకవేళ ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ఒకవేళ (ldc < MAX(1, m))  info = 8;
    } లేకపోతే {
      ఒకవేళ (ldc < MAX(1, n))  info = 8;
    }

    ఒకవేళ ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ఒకవేళ (lda < MAX(1, m))  info = 5;
    } లేకపోతే {
      ఒకవేళ (lda < MAX(1, n))  info = 5;
    }
    ఒకవేళ (n < 0)		  info = 1;
    ఒకవేళ (m < 0)		  info = 2;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  ఒకవేళ ((m==0) || (n==0)) తిరిగి;


  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();


  GEADD_K(m,n,alpha, a, lda, beta, c, ldc,
#నిర్వచితమైతే CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#లేకపోతే
          transa, transc
#అయితే_ముగింపు
          ); 


  FUNCTION_PROFILE_చివర(1, 2* m * n ,  2 * m * n);

  IDEBUG_చివర;

  తిరిగి;

}
