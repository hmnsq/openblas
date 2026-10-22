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

#اگر_تعریف_نشده TRMM
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DTRSM "
#وإلا
#تعریف_کن ERROR_NAME "STRSM "
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZTRSM "
#وإلا
#تعریف_کن ERROR_NAME "CTRSM "
#پایان_اگر
#پایان_اگر
#وإلا
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DTRMM "
#وإلا
#تعریف_کن ERROR_NAME "STRMM "
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZTRMM "
#وإلا
#تعریف_کن ERROR_NAME "CTRMM "
#پایان_اگر
#پایان_اگر
#پایان_اگر

#اگر_تعریف_نشده COMPLEX
#تعریف_کن SMP_ضریب 256
#وإلا
#تعریف_کن SMP_ضریب 128
#پایان_اگر

#اگر_تعریف_شده DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#پایان_اگر                                                               

ساكن صحيح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_تعریف_نشده TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#وإلا
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#پایان_اگر
};

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

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
  صحيح side;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  آرگومانها.m = *M;
  آرگومانها.n = *N;

  آرگومانها.a = (فراغ *)a;
  آرگومانها.b = (فراغ *)b;

  آرگومانها.lda = *ldA;
  آرگومانها.ldb = *ldB;

  آرگومانها.beta = (فراغ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (diag_arg  == 'U') unit  = 0;
  إذا (diag_arg  == 'N') unit  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  nrowa = آرگومانها.m;
  إذا (side & 1) nrowa = آرگومانها.n;

  info = 0;

  إذا (آرگومانها.ldb < MAX(1,آرگومانها.m)) info = 11;
  إذا (آرگومانها.lda < MAX(1,nrowa))  info =  9;
  إذا (آرگومانها.n < 0)               info =  6;
  إذا (آرگومانها.m < 0)               info =  5;
  إذا (unit < 0)                 info =  4;
  إذا (trans < 0)                info =  3;
  إذا (uplo  < 0)                info =  2;
  إذا (side  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME)-1);
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   تعداد CBLAS_SIDE Side,  تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE Trans, تعداد CBLAS_DIAG Diag,
	   blasint m, blasint n,
#اگر_تعریف_نشده COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#پایان_اگر

  blas_arg_t آرگومانها;
  صحيح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *میانگیر;
  XFLOAT *sa, *sb;

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

  آرگومانها.a = (فراغ *)a;
  آرگومانها.b = (فراغ *)b;

  آرگومانها.lda = lda;
  آرگومانها.ldb = ldb;

#اگر_تعریف_نشده COMPLEX
  آرگومانها.beta = (فراغ *)&alpha;
#وإلا
  آرگومانها.beta = (فراغ *)alpha;
#پایان_اگر

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    آرگومانها.m = m;
    آرگومانها.n = n;

    إذا (Side == Cblas_چپ)         side  = 0;
    إذا (Side == Cblas_راست)        side  = 1;

    إذا (Uplo == Cblas_فوقانی)        uplo  = 0;
    إذا (Uplo == Cblas_پایینی)        uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#اگر_تعریف_نشده COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#پایان_اگر

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = آرگومانها.m;
    إذا (side & 1) nrowa = آرگومانها.n;

    إذا (آرگومانها.ldb < MAX(1,آرگومانها.m)) info = 11;
    إذا (آرگومانها.lda < MAX(1,nrowa))  info =  9;
    إذا (آرگومانها.n < 0)               info =  6;
    إذا (آرگومانها.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {
    آرگومانها.m = n;
    آرگومانها.n = m;

    إذا (Side == Cblas_چپ)         side  = 1;
    إذا (Side == Cblas_راست)        side  = 0;

    إذا (Uplo == Cblas_فوقانی)        uplo  = 1;
    إذا (Uplo == Cblas_پایینی)        uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#اگر_تعریف_نشده COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#پایان_اگر

    إذا (Diag == CblasUnit)         unit  = 0;
    إذا (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = آرگومانها.m;
    إذا (side & 1) nrowa = آرگومانها.n;

    إذا (آرگومانها.ldb < MAX(1,آرگومانها.m)) info = 11;
    إذا (آرگومانها.lda < MAX(1,nrowa))  info =  9;
    إذا (آرگومانها.n < 0)               info =  6;
    إذا (آرگومانها.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (آرگومانها.m == 0 || آرگومانها.n == 0) أرجع;

#إذا defined(TRMM)
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#پایان_اگر
)
#پایان_اگر
  إذا (ترتیب == Cblas_ردیف_Major && Diag == CblasNonUnit && Side == Cblas_چپ && m == lda && n == ldb) {
    إذا (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_فوقانی ? STRMM_مستقیم_LNUN : STRMM_مستقیم_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا إذا (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_فوقانی ? STRMM_مستقیم_LTUN : STRMM_مستقیم_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#پایان_اگر
#پایان_اگر
#پایان_اگر



//end of the ifndef CBLAS ... else ... section	
#پایان_اگر

  إذا ((آرگومانها.m == 0) || (آرگومانها.n == 0)) أرجع;
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
  نما |= (trans << BLAS_TRANSA_لغزش);
  نما |= (side  << BLAS_RSIDE_لغزش);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  إذا ( آرگومانها.m * آرگومانها.n < SMP_ضریب * GEMM_MULTITHREAD_آستانه)
	آرگومانها.nthreads = 1;
  وإلا
	آرگومانها.nthreads = رقم_cpu_avail(3);
		

  إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&آرگومانها, عدم, عدم, sa, sb, 0);

#اگر_تعریف_شده SMP
  } وإلا {
    إذا (!side) {
      gemm_thread_n(نما, &آرگومانها, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, آرگومانها.nthreads);
    } وإلا {
      gemm_thread_m(نما, &آرگومانها, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, آرگومانها.nthreads);
    }
  }
#پایان_اگر

  blas_حافظه_free(میانگیر);

  FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE,
		       (!side) ? آرگومانها.m * (آرگومانها.m + آرگومانها.n) : آرگومانها.n * (آرگومانها.m + آرگومانها.n),
		       (!side) ? آرگومانها.m * آرگومانها.m * آرگومانها.n : آرگومانها.m * آرگومانها.n * آرگومانها.n);

  IDEBUG_انتها;

  أرجع;
}

