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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QSYR2K"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DSYR2K"
#وإلا
#وصف_ڪريو ERROR_NAME "SSYR2K"
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل HEMM
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XSYR2K"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZSYR2K"
#وإلا
#وصف_ڪريو ERROR_NAME "CSYR2K"
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XHER2K"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHER2K"
#وإلا
#وصف_ڪريو ERROR_NAME "CHER2K"
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

ساكن صحيح (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جيڪڏهن_وصف_نه_ٿيل HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#وإلا
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#جيڪڏهن_پورو
};

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t دليلون;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دليلون.n = *N;
  دليلون.k = *K;

  دليلون.a = (فراغ *)a;
  دليلون.b = (فراغ *)b;
  دليلون.c = (فراغ *)c;

  دليلون.lda = *ldA;
  دليلون.ldb = *ldB;
  دليلون.ldc = *ldC;

  دليلون.alpha = (فراغ *)alpha;
  دليلون.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
#جيڪڏهن_وصف_ٿيل HEMM
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
#جيڪڏهن_پورو

#جيڪڏهن_پورو


  nrowa = دليلون.n;
  إذا (trans & 1) nrowa = دليلون.k;

  info = 0;

  إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 12;
  إذا (دليلون.ldb < MAX(1,nrowa))  info =  9;
  إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
  إذا (دليلون.k < 0)               info =  4;
  إذا (دليلون.n < 0)               info =  3;
  إذا (trans < 0)                info =  2;
  إذا (uplo  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
#جيڪڏهن_پورو
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#جيڪڏهن_پورو
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	   FLOAT *c,
#وإلا
	   فراغ *vc,
#جيڪڏهن_پورو
           blasint ldc) {

#جيڪڏهن_وصف_ٿيل COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#إذا !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#جيڪڏهن_پورو
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#جيڪڏهن_پورو

  blas_arg_t دليلون;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل HEMM
  FLOAT CAlpha[2];
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  دليلون.n = n;
  دليلون.k = k;

  دليلون.a = (فراغ *)a;
  دليلون.b = (فراغ *)b;
  دليلون.c = (فراغ *)c;

  دليلون.lda = lda;
  دليلون.ldb = ldb;
  دليلون.ldc = ldc;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  دليلون.alpha = (فراغ *)&alpha;
#وإلا
  دليلون.alpha = (فراغ *)alpha;
#جيڪڏهن_پورو

#إذا !defined(COMPLEX) || defined(HEMM)
  دليلون.beta  = (فراغ *)&beta;
#وإلا
  دليلون.beta  = (فراغ *)beta;
#جيڪڏهن_پورو

  trans = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (Uplo == Cblas_مٿانهون) uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasTrans)       trans = 1;
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 1;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 1;
#جيڪڏهن_پورو

    info = -1;

    nrowa = دليلون.n;
    إذا (trans & 1) nrowa = دليلون.k;

    إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 12;
    إذا (دليلون.ldb < MAX(1,nrowa))  info =  9;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
    إذا (دليلون.k < 0)               info =  4;
    إذا (دليلون.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {

#جيڪڏهن_وصف_ٿيل HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    دليلون.alpha = (فراغ *)CAlpha;
#جيڪڏهن_پورو

    إذا (Uplo == Cblas_مٿانهون) uplo  = 1;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasTrans)       trans = 0;
    إذا (Trans == CblasConjNoTrans) trans = 1;
    إذا (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 0;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 0;
#جيڪڏهن_پورو

    info = -1;

    nrowa = دليلون.n;
    إذا (trans & 1) nrowa = دليلون.k;

    إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 12;
    إذا (دليلون.ldb < MAX(1,nrowa))  info =  9;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
    إذا (دليلون.k < 0)               info =  4;
    إذا (دليلون.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دليلون.n == 0) أرجع;

#جيڪڏهن_وصف_ٿيل DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#جيڪڏهن_پورو

#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYR2_K_KERNEL_سڌو)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جيڪڏهن_پورو
)
#جيڪڏهن_پورو
   إذا (ترتيب == Cblas_صف_Major && n == ldc) {
     إذا (Trans == CblasNoTrans && k == lda && k == ldb) {
      إذا (Uplo == Cblas_مٿانهون) {
        SSYR2_K_سڌو_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_هيٺاهون) {
        SSYR2_K_سڌو_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
     وإلا إذا (Trans == CblasTrans && n == lda && n ==ldb) {
      إذا (Uplo == Cblas_مٿانهون) {
        SSYR2_K_سڌو_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_هيٺاهون) {
        SSYR2_K_سڌو_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
   }
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_پورو


  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  بفر = (FLOAT *)blas_يادگيري_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);

#جيڪڏهن_وصف_ٿيل SMP
  إذا (!trans){
    گهڻائي |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    گهڻائي |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  گهڻائي |= (uplo  << BLAS_UPLO_سرڪايو);

  دليلون.common = عدم;
  إذا (دليلون.n*دليلون.k <1000)
	  دليلون.nthreads =1 ;
  وإلا
  دليلون.nthreads = انگ_cpu_avail(3);

  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    (syr2k[(uplo << 1) | trans ])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP

  } وإلا {

    syrk_thread(گهڻائي, &دليلون, عدم, عدم, syr2k[(uplo << 1) | trans ], sa, sb, دليلون.nthreads);

  }
#جيڪڏهن_پورو

  blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE, 2 * دليلون.n * دليلون.k + دليلون.n * دليلون.n, 2 * دليلون.n * دليلون.n * دليلون.k);

  IDEBUG_پڇاڙي;

  أرجع;
}
