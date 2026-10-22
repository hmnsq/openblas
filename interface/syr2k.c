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
#ضمن <ctype.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QSYR2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DSYR2K"
#وإلا
#تعریف_کن ERROR_NAME "SSYR2K"
#پایان_اگر
#وإلا
#اگر_تعریف_نشده HEMM
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XSYR2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZSYR2K"
#وإلا
#تعریف_کن ERROR_NAME "CSYR2K"
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XHER2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZHER2K"
#وإلا
#تعریف_کن ERROR_NAME "CHER2K"
#پایان_اگر
#پایان_اگر
#پایان_اگر

ساكن صحيح (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_تعریف_نشده HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#وإلا
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#پایان_اگر
};

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t آرگومانها;

  FLOAT *میانگیر;
  FLOAT *sa, *sb;

#اگر_تعریف_شده SMP
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
  صحيح نما  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح نما  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح نما  =  BLAS_SINGLE  | BLAS_REAL;
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
  صحيح نما  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح نما  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح نما  =  BLAS_SINGLE  | BLAS_COMPLEX;
#پایان_اگر
#پایان_اگر
#پایان_اگر

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  آرگومانها.n = *N;
  آرگومانها.k = *K;

  آرگومانها.a = (فراغ *)a;
  آرگومانها.b = (فراغ *)b;
  آرگومانها.c = (فراغ *)c;

  آرگومانها.lda = *ldA;
  آرگومانها.ldb = *ldB;
  آرگومانها.ldc = *ldC;

  آرگومانها.alpha = (فراغ *)alpha;
  آرگومانها.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

#اگر_تعریف_نشده COMPLEX
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
#اگر_تعریف_شده HEMM
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
#پایان_اگر

#پایان_اگر


  nrowa = آرگومانها.n;
  إذا (trans & 1) nrowa = آرگومانها.k;

  info = 0;

  إذا (آرگومانها.ldc < MAX(1,آرگومانها.n)) info = 12;
  إذا (آرگومانها.ldb < MAX(1,nrowa))  info =  9;
  إذا (آرگومانها.lda < MAX(1,nrowa))  info =  7;
  إذا (آرگومانها.k < 0)               info =  4;
  إذا (آرگومانها.n < 0)               info =  3;
  إذا (trans < 0)                info =  2;
  إذا (uplo  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#اگر_تعریف_نشده COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
#پایان_اگر
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#پایان_اگر
#اگر_تعریف_نشده COMPLEX
	   FLOAT *c,
#وإلا
	   فراغ *vc,
#پایان_اگر
           blasint ldc) {

#اگر_تعریف_شده COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#إذا !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#پایان_اگر
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#پایان_اگر

  blas_arg_t آرگومانها;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *میانگیر;
  FLOAT *sa, *sb;

#اگر_تعریف_شده HEMM
  FLOAT CAlpha[2];
#پایان_اگر

#اگر_تعریف_شده SMP
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
  صحيح نما  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح نما  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح نما  =  BLAS_SINGLE  | BLAS_REAL;
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
  صحيح نما  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح نما  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح نما  =  BLAS_SINGLE  | BLAS_COMPLEX;
#پایان_اگر
#پایان_اگر
#پایان_اگر

  PRINT_DEBUG_CNAME;

  آرگومانها.n = n;
  آرگومانها.k = k;

  آرگومانها.a = (فراغ *)a;
  آرگومانها.b = (فراغ *)b;
  آرگومانها.c = (فراغ *)c;

  آرگومانها.lda = lda;
  آرگومانها.ldb = ldb;
  آرگومانها.ldc = ldc;

#اگر_تعریف_نشده COMPLEX
  آرگومانها.alpha = (فراغ *)&alpha;
#وإلا
  آرگومانها.alpha = (فراغ *)alpha;
#پایان_اگر

#إذا !defined(COMPLEX) || defined(HEMM)
  آرگومانها.beta  = (فراغ *)&beta;
#وإلا
  آرگومانها.beta  = (فراغ *)beta;
#پایان_اگر

  trans = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_فوقانی) uplo  = 0;
    إذا (Uplo == Cblas_پایینی) uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
#اگر_تعریف_نشده COMPLEX
    إذا (Trans == CblasTrans)       trans = 1;
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 1;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 1;
#پایان_اگر

    info = -1;

    nrowa = آرگومانها.n;
    إذا (trans & 1) nrowa = آرگومانها.k;

    إذا (آرگومانها.ldc < MAX(1,آرگومانها.n)) info = 12;
    إذا (آرگومانها.ldb < MAX(1,nrowa))  info =  9;
    إذا (آرگومانها.lda < MAX(1,nrowa))  info =  7;
    إذا (آرگومانها.k < 0)               info =  4;
    إذا (آرگومانها.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {

#اگر_تعریف_شده HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    آرگومانها.alpha = (فراغ *)CAlpha;
#پایان_اگر

    إذا (Uplo == Cblas_فوقانی) uplo  = 1;
    إذا (Uplo == Cblas_پایینی) uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 1;
#اگر_تعریف_نشده COMPLEX
    إذا (Trans == CblasTrans)       trans = 0;
    إذا (Trans == CblasConjNoTrans) trans = 1;
    إذا (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 0;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 0;
#پایان_اگر

    info = -1;

    nrowa = آرگومانها.n;
    إذا (trans & 1) nrowa = آرگومانها.k;

    إذا (آرگومانها.ldc < MAX(1,آرگومانها.n)) info = 12;
    إذا (آرگومانها.ldb < MAX(1,nrowa))  info =  9;
    إذا (آرگومانها.lda < MAX(1,nrowa))  info =  7;
    إذا (آرگومانها.k < 0)               info =  4;
    إذا (آرگومانها.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (آرگومانها.n == 0) أرجع;

#اگر_تعریف_شده DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#پایان_اگر

#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYR2_K_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#پایان_اگر
)
#پایان_اگر
   إذا (ترتیب == Cblas_ردیف_Major && n == ldc) {
     إذا (Trans == CblasNoTrans && k == lda && k == ldb) {
      إذا (Uplo == Cblas_فوقانی) {
        SSYR2_K_مستقیم_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_پایینی) {
        SSYR2_K_مستقیم_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
     وإلا إذا (Trans == CblasTrans && n == lda && n ==ldb) {
      إذا (Uplo == Cblas_فوقانی) {
        SSYR2_K_مستقیم_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_پایینی) {
        SSYR2_K_مستقیم_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
   }
#پایان_اگر
#پایان_اگر

#پایان_اگر


  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  میانگیر = (FLOAT *)blas_حافظه_alloc(0);
  إذا (!میانگیر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);

#اگر_تعریف_شده SMP
  إذا (!trans){
    نما |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    نما |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  نما |= (uplo  << BLAS_UPLO_لغزش);

  آرگومانها.common = عدم;
  إذا (آرگومانها.n*آرگومانها.k <1000)
	  آرگومانها.nthreads =1 ;
  وإلا
  آرگومانها.nthreads = رقم_cpu_avail(3);

  إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

    (syr2k[(uplo << 1) | trans ])(&آرگومانها, عدم, عدم, sa, sb, 0);

#اگر_تعریف_شده SMP

  } وإلا {

    syrk_thread(نما, &آرگومانها, عدم, عدم, syr2k[(uplo << 1) | trans ], sa, sb, آرگومانها.nthreads);

  }
#پایان_اگر

  blas_حافظه_free(میانگیر);

  FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE, 2 * آرگومانها.n * آرگومانها.k + آرگومانها.n * آرگومانها.n, 2 * آرگومانها.n * آرگومانها.n * آرگومانها.k);

  IDEBUG_انتها;

  أرجع;
}
