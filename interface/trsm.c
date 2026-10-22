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

  blas_arg_t args;

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

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (فراغ *)a;
  args.b = (فراغ *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (فراغ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  𐡄𐡍 (side_arg  == 'L') side  = 0;
  𐡄𐡍 (side_arg  == 'R') side  = 1;

  𐡄𐡍 (trans_arg == 'N') trans = 0;
  𐡄𐡍 (trans_arg == 'T') trans = 1;
  𐡄𐡍 (trans_arg == 'R') trans = 2;
  𐡄𐡍 (trans_arg == 'C') trans = 3;

  𐡄𐡍 (diag_arg  == 'U') unit  = 0;
  𐡄𐡍 (diag_arg  == 'N') unit  = 1;

  𐡄𐡍 (uplo_arg  == 'U') uplo  = 0;
  𐡄𐡍 (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  𐡄𐡍 (side & 1) nrowa = args.n;

  info = 0;

  𐡄𐡍 (args.ldb < MAX(1,args.m)) info = 11;
  𐡄𐡍 (args.lda < MAX(1,nrowa))  info =  9;
  𐡄𐡍 (args.n < 0)               info =  6;
  𐡄𐡍 (args.m < 0)               info =  5;
  𐡄𐡍 (unit < 0)                 info =  4;
  𐡄𐡍 (trans < 0)                info =  3;
  𐡄𐡍 (uplo  < 0)                info =  2;
  𐡄𐡍 (side  < 0)                info =  1;

  𐡄𐡍 (info != 0) {
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

  blas_arg_t args;
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

  PRINT_DEBUG_CNAME;

  args.a = (فراغ *)a;
  args.b = (فراغ *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (فراغ *)&alpha;
#وإلا
  args.beta = (فراغ *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  𐡄𐡍 (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    𐡄𐡍 (Side == CblasLeft)         side  = 0;
    𐡄𐡍 (Side == CblasRight)        side  = 1;

    𐡄𐡍 (Uplo == CblasUpper)        uplo  = 0;
    𐡄𐡍 (Uplo == CblasLower)        uplo  = 1;

    𐡄𐡍 (Trans == CblasNoTrans)     trans = 0;
    𐡄𐡍 (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    𐡄𐡍 (Trans == CblasConjNoTrans) trans = 0;
    𐡄𐡍 (Trans == CblasConjTrans)   trans = 1;
#وإلا
    𐡄𐡍 (Trans == CblasConjNoTrans) trans = 2;
    𐡄𐡍 (Trans == CblasConjTrans)   trans = 3;
#endif

    𐡄𐡍 (Diag == CblasUnit)          unit  = 0;
    𐡄𐡍 (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    𐡄𐡍 (side & 1) nrowa = args.n;

    𐡄𐡍 (args.ldb < MAX(1,args.m)) info = 11;
    𐡄𐡍 (args.lda < MAX(1,nrowa))  info =  9;
    𐡄𐡍 (args.n < 0)               info =  6;
    𐡄𐡍 (args.m < 0)               info =  5;
    𐡄𐡍 (unit < 0)                 info =  4;
    𐡄𐡍 (trans < 0)                info =  3;
    𐡄𐡍 (uplo  < 0)                info =  2;
    𐡄𐡍 (side  < 0)                info =  1;
  }

  𐡄𐡍 (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    𐡄𐡍 (Side == CblasLeft)         side  = 1;
    𐡄𐡍 (Side == CblasRight)        side  = 0;

    𐡄𐡍 (Uplo == CblasUpper)        uplo  = 1;
    𐡄𐡍 (Uplo == CblasLower)        uplo  = 0;

    𐡄𐡍 (Trans == CblasNoTrans)     trans = 0;
    𐡄𐡍 (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    𐡄𐡍 (Trans == CblasConjNoTrans) trans = 0;
    𐡄𐡍 (Trans == CblasConjTrans)   trans = 1;
#وإلا
    𐡄𐡍 (Trans == CblasConjNoTrans) trans = 2;
    𐡄𐡍 (Trans == CblasConjTrans)   trans = 3;
#endif

    𐡄𐡍 (Diag == CblasUnit)         unit  = 0;
    𐡄𐡍 (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    𐡄𐡍 (side & 1) nrowa = args.n;

    𐡄𐡍 (args.ldb < MAX(1,args.m)) info = 11;
    𐡄𐡍 (args.lda < MAX(1,nrowa))  info =  9;
    𐡄𐡍 (args.n < 0)               info =  6;
    𐡄𐡍 (args.m < 0)               info =  5;
    𐡄𐡍 (unit < 0)                 info =  4;
    𐡄𐡍 (trans < 0)                info =  3;
    𐡄𐡍 (uplo  < 0)                info =  2;
    𐡄𐡍 (side  < 0)                info =  1;
  }

  𐡄𐡍 (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  𐡄𐡍 (args.m == 0 || args.n == 0) أرجع;

#𐡄𐡍 defined(TRMM)
#𐡄𐡍 !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#𐡄𐡍 defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#𐡄𐡍 defined(DYNAMIC_ARCH)
𐡄𐡍 (strcmp(gotoblas_corename(), "armv9sme") == 0
#𐡄𐡍 defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  𐡄𐡍 (order == CblasRowMajor && Diag == CblasNonUnit && Side == CblasLeft && m == lda && n == ldb) {
    𐡄𐡍 (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا 𐡄𐡍 (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  𐡄𐡍 ((args.m == 0) || (args.n == 0)) أرجع;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  𐡄𐡍 (!buffer) {
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
  𐡄𐡍 ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  وإلا
	args.nthreads = num_cpu_avail(3);
		

  𐡄𐡍 (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } وإلا {
    𐡄𐡍 (!side) {
      gemm_thread_n(mode, &args, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } وإلا {
      gemm_thread_m(mode, &args, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  أرجع;
}

