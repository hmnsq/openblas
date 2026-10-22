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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#وإلا
#define ERROR_NAME "STRSM "
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#وإلا
#define ERROR_NAME "CTRSM "
#endif
#endif
#وإلا
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#وإلا
#define ERROR_NAME "STRMM "
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#وإلا
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#وإلا
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#endif                                                               

ساكن صحيح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
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
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blas_arg_t حجوجي;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وإلا
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  صحيح side;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  صحيح nrowa;

  بغ_DEBUG_NAME;

  حجوجي.m = *M;
  حجوجي.n = *N;

  حجوجي.a = (فراغ *)a;
  حجوجي.b = (فراغ *)b;

  حجوجي.lda = *ldA;
  حجوجي.ldb = *ldB;

  حجوجي.beta = (فراغ *)alpha;

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

  nrowa = حجوجي.m;
  إذا (side & 1) nrowa = حجوجي.n;

  info = 0;

  إذا (حجوجي.ldb < MAX(1,حجوجي.m)) info = 11;
  إذا (حجوجي.lda < MAX(1,nrowa))  info =  9;
  إذا (حجوجي.n < 0)               info =  6;
  إذا (حجوجي.m < 0)               info =  5;
  إذا (unit < 0)                 info =  4;
  إذا (trans < 0)                info =  3;
  إذا (uplo  < 0)                info =  2;
  إذا (side  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME)-1);
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order,
	   تعداد CBLAS_SIDE Side,  تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE Trans, تعداد CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
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
#endif

  blas_arg_t حجوجي;
  صحيح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وإلا
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  بغ_DEBUG_CNAME;

  حجوجي.a = (فراغ *)a;
  حجوجي.b = (فراغ *)b;

  حجوجي.lda = lda;
  حجوجي.ldb = ldb;

#ifndef COMPLEX
  حجوجي.beta = (فراغ *)&alpha;
#وإلا
  حجوجي.beta = (فراغ *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  إذا (order == CblasColMajor) {
    حجوجي.m = m;
    حجوجي.n = n;

    إذا (Side == Cblas_هغو)         side  = 0;
    إذا (Side == Cblas_دام)        side  = 1;

    إذا (Uplo == CblasUpper)        uplo  = 0;
    إذا (Uplo == CblasLower)        uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#endif

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = حجوجي.m;
    إذا (side & 1) nrowa = حجوجي.n;

    إذا (حجوجي.ldb < MAX(1,حجوجي.m)) info = 11;
    إذا (حجوجي.lda < MAX(1,nrowa))  info =  9;
    إذا (حجوجي.n < 0)               info =  6;
    إذا (حجوجي.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (order == CblasRowMajor) {
    حجوجي.m = n;
    حجوجي.n = m;

    إذا (Side == Cblas_هغو)         side  = 1;
    إذا (Side == Cblas_دام)        side  = 0;

    إذا (Uplo == CblasUpper)        uplo  = 1;
    إذا (Uplo == CblasLower)        uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#endif

    إذا (Diag == CblasUnit)         unit  = 0;
    إذا (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = حجوجي.m;
    إذا (side & 1) nrowa = حجوجي.n;

    إذا (حجوجي.ldb < MAX(1,حجوجي.m)) info = 11;
    إذا (حجوجي.lda < MAX(1,nrowa))  info =  9;
    إذا (حجوجي.n < 0)               info =  6;
    إذا (حجوجي.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (حجوجي.m == 0 || حجوجي.n == 0) أرجع;

#إذا defined(TRMM)
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  إذا (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_هغو && m == lda && n == ldb) {
    إذا (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا إذا (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  إذا ((حجوجي.m == 0) || (حجوجي.n == 0)) أرجع;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  إذا (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  mode |= (trans << BLAS_TRANSA_SHIFT);
  mode |= (side  << BLAS_RSIDE_SHIFT);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  إذا ( حجوجي.m * حجوجي.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	حجوجي.nthreads = 1;
  وإلا
	حجوجي.nthreads = لمب_cpu_avail(3);
		

  إذا (حجوجي.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&حجوجي, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } وإلا {
    إذا (!side) {
      gemm_thread_n(mode, &حجوجي, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, حجوجي.nthreads);
    } وإلا {
      gemm_thread_m(mode, &حجوجي, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, حجوجي.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? حجوجي.m * (حجوجي.m + حجوجي.n) : حجوجي.n * (حجوجي.m + حجوجي.n),
		       (!side) ? حجوجي.m * حجوجي.m * حجوجي.n : حجوجي.m * حجوجي.n * حجوجي.n);

  IDEBUG_END;

  أرجع;
}

