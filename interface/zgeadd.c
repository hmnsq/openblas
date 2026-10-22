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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEADD "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEADD "
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(blasint *M, blasint *N, FLOAT *ALPHA, FLOAT *a, blasint *LDA,
		                  FLOAT *BETA,  FLOAT *c, blasint *LDC,ಅಕ್ಷರ *TRANS_A,ಅಕ್ಷರ *TRANS_C)
{

  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint ldc = *LDC; 

  blasint info;
ಅಕ್ಷರ transa = (*TRANS_A == 'T') || (*TRANS_A == 't') || (*TRANS_A == 'C') || (*TRANS_A == 'c');
ಅಕ್ಷರ transc = (*TRANS_C == 'T') || (*TRANS_C == 't') || (*TRANS_C == 'C') || (*TRANS_C == 'c');


  PRINT_DEBUG_NAME;

  info = 0;


ಆದರೆ(transa){
  ಆದರೆ (lda < MAX(1, n))info = 5;


  }	ಇಲ್ಲದಿದ್ದರೆ{

   ಆದರೆ (lda < MAX(1, m))	info = 5;
  }
 

  ಆದರೆ(transc){
        ಆದರೆ (ldc < MAX(1, n))	info = 8;


   }ಇಲ್ಲದಿದ್ದರೆ{
       
      ಆದರೆ (ldc < MAX(1, m))	info = 8;
          }


  ಆದರೆ (n < 0)		info = 2;
  ಆದರೆ (m < 0)		info = 1;

  ಆದರೆ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ
ಶೂನ್ಯ CNAME( ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,ಎಣಿಕೆ CBLAS_TRANSPOSE transa, ಎಣಿಕೆ CBLAS_TRANSPOSE transc, blasint m,  blasint n,  FLOAT *ALPHA, FLOAT *a,  blasint lda, FLOAT *BETA, 
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

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {

    info = -1;

ಆದರೆ ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ಆದರೆ (ldc < MAX(1, m))  info = 8;
    } ಇಲ್ಲದಿದ್ದರೆ {
      ಆದರೆ (ldc < MAX(1, n))  info = 8;
    }

    ಆದರೆ ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ಆದರೆ (lda < MAX(1, m))  info = 5;
    } ಇಲ್ಲದಿದ್ದರೆ {
      ಆದರೆ (lda < MAX(1, n))  info = 5;
    }
    ಆದರೆ (n < 0)		  info = 2;
    ಆದರೆ (m < 0)		  info = 1;  
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    info = -1;

    t = n;
    n = m;
    m = t;
ಆದರೆ ( (transc == CblasNoTrans) || (transc == CblasConjNoTrans) ) {
      ಆದರೆ (ldc < MAX(1, m))  info = 8;
    } ಇಲ್ಲದಿದ್ದರೆ {
      ಆದರೆ (ldc < MAX(1, n))  info = 8;
    }

    ಆದರೆ ( (transa == CblasNoTrans) || (transa == CblasConjNoTrans) ) {
      ಆದರೆ (lda < MAX(1, m))  info = 5;
    } ಇಲ್ಲದಿದ್ದರೆ {
      ಆದರೆ (lda < MAX(1, n))  info = 5;
    }
    ಆದರೆ (n < 0)		  info = 1;
    ಆದರೆ (m < 0)		  info = 2;
  }

 
  

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ ((m==0) || (n==0)) ಹಿಂತಿರುಗಿ;


  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();


  GEADD_K(m,n,ALPHA[0],ALPHA[1], a, lda, BETA[0], BETA[1], c, ldc,
        
          #ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ CBLAS
          (transa == CblasTrans) || (transa == CblasConjTrans),
          (transc == CblasTrans) || (transc == CblasConjTrans)
#ಇಲ್ಲದಿದ್ದರೆ
          transa, transc
#ಆದರೆ_ಅಂತ್ಯ
          );


           


  FUNCTION_PROFILE_ಕೊನೆ(1, 2* m * n ,  2 * m * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;

}
