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
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QSYR2K"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DSYR2K"
#अन्यथा
#परिभाषय ERROR_NAME "SSYR2K"
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् HEMM
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XSYR2K"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZSYR2K"
#अन्यथा
#परिभाषय ERROR_NAME "CSYR2K"
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XHER2K"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZHER2K"
#अन्यथा
#परिभाषय ERROR_NAME "CHER2K"
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

स्थिर पूर्णाङ्क (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#यदि_अपरिभाषितम् HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#अन्यथा
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#यदि_अन्तः
};

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;

  blas_arg_t तर्काः;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क trans;
  पूर्णाङ्क nrowa;

  PRINT_DEBUG_NAME;

  तर्काः.n = *N;
  तर्काः.k = *K;

  तर्काः.a = (शून्य *)a;
  तर्काः.b = (शून्य *)b;
  तर्काः.c = (शून्य *)c;

  तर्काः.lda = *ldA;
  तर्काः.ldb = *ldB;
  तर्काः.ldc = *ldC;

  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

#यदि_अपरिभाषितम् COMPLEX
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
#यदि_परिभाषितम् HEMM
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
#यदि_अन्तः

#यदि_अन्तः


  nrowa = तर्काः.n;
  यदि (trans & 1) nrowa = तर्काः.k;

  info = 0;

  यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 12;
  यदि (तर्काः.ldb < MAX(1,nrowa))  info =  9;
  यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
  यदि (तर्काः.k < 0)               info =  4;
  यदि (तर्काः.n < 0)               info =  3;
  यदि (trans < 0)                info =  2;
  यदि (uplo  < 0)                info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo, गणना CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
#यदि_अन्तः
#यदि !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   शून्य *vbeta,
#यदि_अन्तः
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT *c,
#अन्यथा
	   शून्य *vc,
#यदि_अन्तः
           blasint ldc) {

#यदि_परिभाषितम् COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#यदि !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#यदि_अन्तः
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#यदि_अन्तः

  blas_arg_t तर्काः;
  पूर्णाङ्क uplo, trans;
  blasint info, nrowa;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि_परिभाषितम् HEMM
  FLOAT CAlpha[2];
#यदि_अन्तः

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  तर्काः.n = n;
  तर्काः.k = k;

  तर्काः.a = (शून्य *)a;
  तर्काः.b = (शून्य *)b;
  तर्काः.c = (शून्य *)c;

  तर्काः.lda = lda;
  तर्काः.ldb = ldb;
  तर्काः.ldc = ldc;

#यदि_अपरिभाषितम् COMPLEX
  तर्काः.alpha = (शून्य *)&alpha;
#अन्यथा
  तर्काः.alpha = (शून्य *)alpha;
#यदि_अन्तः

#यदि !defined(COMPLEX) || defined(HEMM)
  तर्काः.beta  = (शून्य *)&beta;
#अन्यथा
  तर्काः.beta  = (शून्य *)beta;
#यदि_अन्तः

  trans = -1;
  uplo  = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 1;

    यदि (Trans == CblasNoTrans)     trans = 0;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasTrans)       trans = 1;
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 1;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 1;
#यदि_अन्तः

    info = -1;

    nrowa = तर्काः.n;
    यदि (trans & 1) nrowa = तर्काः.k;

    यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 12;
    यदि (तर्काः.ldb < MAX(1,nrowa))  info =  9;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
    यदि (तर्काः.k < 0)               info =  4;
    यदि (तर्काः.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {

#यदि_परिभाषितम् HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    तर्काः.alpha = (शून्य *)CAlpha;
#यदि_अन्तः

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 0;

    यदि (Trans == CblasNoTrans)     trans = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasTrans)       trans = 0;
    यदि (Trans == CblasConjNoTrans) trans = 1;
    यदि (Trans == CblasConjTrans)   trans = 0;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 0;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 0;
#यदि_अन्तः

    info = -1;

    nrowa = तर्काः.n;
    यदि (trans & 1) nrowa = तर्काः.k;

    यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 12;
    यदि (तर्काः.ldb < MAX(1,nrowa))  info =  9;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
    यदि (तर्काः.k < 0)               info =  4;
    यदि (तर्काः.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (तर्काः.n == 0) प्रत्यागम;

#यदि_परिभाषितम् DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#यदि_अन्तः

#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SSYR2_K_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः
   यदि (क्रमः == CblasRowMajor && n == ldc) {
     यदि (Trans == CblasNoTrans && k == lda && k == ldb) {
      यदि (Uplo == Cblas_ऊर्ध्वस्थम्) {
        SSYR2_K_प्रत्यक्षम्_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अन्यथा यदि (Uplo == Cblas_अधस्तनम्) {
        SSYR2_K_प्रत्यक्षम्_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }
     }
     अन्यथा यदि (Trans == CblasTrans && n == lda && n ==ldb) {
      यदि (Uplo == Cblas_ऊर्ध्वस्थम्) {
        SSYR2_K_प्रत्यक्षम्_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अन्यथा यदि (Uplo == Cblas_अधस्तनम्) {
        SSYR2_K_प्रत्यक्षम्_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }
     }
   }
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  बफर = (FLOAT *)blas_स्मृतिः_alloc(0);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);

#यदि_परिभाषितम् SMP
  यदि (!trans){
    बहुलकम् |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    बहुलकम् |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  बहुलकम् |= (uplo  << BLAS_UPLO_सरणम्);

  तर्काः.common = रिक्त;
  यदि (तर्काः.n*तर्काः.k <1000)
	  तर्काः.nthreads =1 ;
  अन्यथा
  तर्काः.nthreads = अङ्कः_cpu_avail(3);

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    (syr2k[(uplo << 1) | trans ])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP

  } अन्यथा {

    syrk_thread(बहुलकम्, &तर्काः, रिक्त, रिक्त, syr2k[(uplo << 1) | trans ], sa, sb, तर्काः.nthreads);

  }
#यदि_अन्तः

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, 2 * तर्काः.n * तर्काः.k + तर्काः.n * तर्काः.n, 2 * तर्काः.n * तर्काः.n * तर्काः.k);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
