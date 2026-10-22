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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QSYR2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DSYR2K"
#وإلا
#تعریف_کړه ERROR_NAME "SSYR2K"
#د_که_پای
#وإلا
#که_نه_تعریف_شوی HEMM
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XSYR2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZSYR2K"
#وإلا
#تعریف_کړه ERROR_NAME "CSYR2K"
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XHER2K"
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZHER2K"
#وإلا
#تعریف_کړه ERROR_NAME "CHER2K"
#د_که_پای
#د_که_پای
#د_که_پای

ساكن صحيح (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#که_نه_تعریف_شوی HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#وإلا
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#د_که_پای
};

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t دلیلونه;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دلیلونه.n = *N;
  دلیلونه.k = *K;

  دلیلونه.a = (فراغ *)a;
  دلیلونه.b = (فراغ *)b;
  دلیلونه.c = (فراغ *)c;

  دلیلونه.lda = *ldA;
  دلیلونه.ldb = *ldB;
  دلیلونه.ldc = *ldC;

  دلیلونه.alpha = (فراغ *)alpha;
  دلیلونه.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

#که_نه_تعریف_شوی COMPLEX
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
#که_تعریف_شوی HEMM
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
#د_که_پای

#د_که_پای


  nrowa = دلیلونه.n;
  إذا (trans & 1) nrowa = دلیلونه.k;

  info = 0;

  إذا (دلیلونه.ldc < MAX(1,دلیلونه.n)) info = 12;
  إذا (دلیلونه.ldb < MAX(1,nrowa))  info =  9;
  إذا (دلیلونه.lda < MAX(1,nrowa))  info =  7;
  إذا (دلیلونه.k < 0)               info =  4;
  إذا (دلیلونه.n < 0)               info =  3;
  إذا (trans < 0)                info =  2;
  إذا (uplo  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#که_نه_تعریف_شوی COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
#د_که_پای
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#د_که_پای
#که_نه_تعریف_شوی COMPLEX
	   FLOAT *c,
#وإلا
	   فراغ *vc,
#د_که_پای
           blasint ldc) {

#که_تعریف_شوی COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#إذا !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#د_که_پای
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#د_که_پای

  blas_arg_t دلیلونه;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#که_تعریف_شوی HEMM
  FLOAT CAlpha[2];
#د_که_پای

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای

  PRINT_DEBUG_CNAME;

  دلیلونه.n = n;
  دلیلونه.k = k;

  دلیلونه.a = (فراغ *)a;
  دلیلونه.b = (فراغ *)b;
  دلیلونه.c = (فراغ *)c;

  دلیلونه.lda = lda;
  دلیلونه.ldb = ldb;
  دلیلونه.ldc = ldc;

#که_نه_تعریف_شوی COMPLEX
  دلیلونه.alpha = (فراغ *)&alpha;
#وإلا
  دلیلونه.alpha = (فراغ *)alpha;
#د_که_پای

#إذا !defined(COMPLEX) || defined(HEMM)
  دلیلونه.beta  = (فراغ *)&beta;
#وإلا
  دلیلونه.beta  = (فراغ *)beta;
#د_که_پای

  trans = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_پاسنی) uplo  = 0;
    إذا (Uplo == Cblas_ښکتنی) uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
#که_نه_تعریف_شوی COMPLEX
    إذا (Trans == CblasTrans)       trans = 1;
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 1;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 1;
#د_که_پای

    info = -1;

    nrowa = دلیلونه.n;
    إذا (trans & 1) nrowa = دلیلونه.k;

    إذا (دلیلونه.ldc < MAX(1,دلیلونه.n)) info = 12;
    إذا (دلیلونه.ldb < MAX(1,nrowa))  info =  9;
    إذا (دلیلونه.lda < MAX(1,nrowa))  info =  7;
    إذا (دلیلونه.k < 0)               info =  4;
    إذا (دلیلونه.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {

#که_تعریف_شوی HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    دلیلونه.alpha = (فراغ *)CAlpha;
#د_که_پای

    إذا (Uplo == Cblas_پاسنی) uplo  = 1;
    إذا (Uplo == Cblas_ښکتنی) uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (Trans == CblasTrans)       trans = 0;
    إذا (Trans == CblasConjNoTrans) trans = 1;
    إذا (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 0;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 0;
#د_که_پای

    info = -1;

    nrowa = دلیلونه.n;
    إذا (trans & 1) nrowa = دلیلونه.k;

    إذا (دلیلونه.ldc < MAX(1,دلیلونه.n)) info = 12;
    إذا (دلیلونه.ldb < MAX(1,nrowa))  info =  9;
    إذا (دلیلونه.lda < MAX(1,nrowa))  info =  7;
    إذا (دلیلونه.k < 0)               info =  4;
    إذا (دلیلونه.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دلیلونه.n == 0) أرجع;

#که_تعریف_شوی DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#د_که_پای

#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYR2_K_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#د_که_پای
)
#د_که_پای
   إذا (ترتیب == Cblas_لیکه_Major && n == ldc) {
     إذا (Trans == CblasNoTrans && k == lda && k == ldb) {
      إذا (Uplo == Cblas_پاسنی) {
        SSYR2_K_مستقیم_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_ښکتنی) {
        SSYR2_K_مستقیم_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
     وإلا إذا (Trans == CblasTrans && n == lda && n ==ldb) {
      إذا (Uplo == Cblas_پاسنی) {
        SSYR2_K_مستقیم_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا إذا (Uplo == Cblas_ښکتنی) {
        SSYR2_K_مستقیم_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
   }
#د_که_پای
#د_که_پای

#د_که_پای


  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  بفر = (FLOAT *)blas_حافظه_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_دننه_واټن_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_دننه_واټن_B);

#که_تعریف_شوی SMP
  إذا (!trans){
    ډېرپېښ |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    ډېرپېښ |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  ډېرپېښ |= (uplo  << BLAS_UPLO_ښويدنه);

  دلیلونه.common = عدم;
  إذا (دلیلونه.n*دلیلونه.k <1000)
	  دلیلونه.nthreads =1 ;
  وإلا
  دلیلونه.nthreads = عدد_cpu_avail(3);

  إذا (دلیلونه.nthreads == 1) {
#د_که_پای

    (syr2k[(uplo << 1) | trans ])(&دلیلونه, عدم, عدم, sa, sb, 0);

#که_تعریف_شوی SMP

  } وإلا {

    syrk_thread(ډېرپېښ, &دلیلونه, عدم, عدم, syr2k[(uplo << 1) | trans ], sa, sb, دلیلونه.nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE * COMPSIZE, 2 * دلیلونه.n * دلیلونه.k + دلیلونه.n * دلیلونه.n, 2 * دلیلونه.n * دلیلونه.n * دلیلونه.k);

  IDEBUG_پای_ټکی;

  أرجع;
}
