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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#অন্যথা
#define ERROR_NAME "STRSM "
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#অন্যথা
#define ERROR_NAME "CTRSM "
#endif
#endif
#অন্যথা
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#অন্যথা
#define ERROR_NAME "STRMM "
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#অন্যথা
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#অন্যথা
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
বাহ্য বর্ণ* gotoblas_corename(শূন্য);                                
#endif                                                               

স্থির পূর্ণাঙ্ক (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#অন্যথা
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

শূন্য NAME(বর্ণ *SIDE, বর্ণ *UPLO, বর্ণ *TRANS, বর্ণ *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  বর্ণ side_arg  = *SIDE;
  বর্ণ uplo_arg  = *UPLO;
  বর্ণ trans_arg = *TRANS;
  বর্ণ diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  পূর্ণাঙ্ক side;
  পূর্ণাঙ্ক uplo;
  পূর্ণাঙ্ক unit;
  পূর্ণাঙ্ক trans;
  পূর্ণাঙ্ক nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (শূন্য *)a;
  args.b = (শূন্য *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (শূন্য *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  যদি (side_arg  == 'L') side  = 0;
  যদি (side_arg  == 'R') side  = 1;

  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'T') trans = 1;
  যদি (trans_arg == 'R') trans = 2;
  যদি (trans_arg == 'C') trans = 3;

  যদি (diag_arg  == 'U') unit  = 0;
  যদি (diag_arg  == 'N') unit  = 1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  যদি (side & 1) nrowa = args.n;

  info = 0;

  যদি (args.ldb < MAX(1,args.m)) info = 11;
  যদি (args.lda < MAX(1,nrowa))  info =  9;
  যদি (args.n < 0)               info =  6;
  যদি (args.m < 0)               info =  5;
  যদি (unit < 0)                 info =  4;
  যদি (trans < 0)                info =  3;
  যদি (uplo  < 0)                info =  2;
  যদি (side  < 0)                info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME)-1);
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_SIDE Side,  গণনা CBLAS_UPLO Uplo,
	   গণনা CBLAS_TRANSPOSE Trans, গণনা CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#অন্যথা
	   শূন্য *valpha,
	   শূন্য *va, blasint lda,
	   শূন্য *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#endif

  blas_arg_t args;
  পূর্ণাঙ্ক side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (শূন্য *)a;
  args.b = (শূন্য *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (শূন্য *)&alpha;
#অন্যথা
  args.beta = (শূন্য *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  যদি (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    যদি (Side == Cblas_বাম)         side  = 0;
    যদি (Side == Cblas_ডান)        side  = 1;

    যদি (Uplo == CblasUpper)        uplo  = 0;
    যদি (Uplo == CblasLower)        uplo  = 1;

    যদি (Trans == CblasNoTrans)     trans = 0;
    যদি (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    যদি (Trans == CblasConjNoTrans) trans = 0;
    যদি (Trans == CblasConjTrans)   trans = 1;
#অন্যথা
    যদি (Trans == CblasConjNoTrans) trans = 2;
    যদি (Trans == CblasConjTrans)   trans = 3;
#endif

    যদি (Diag == CblasUnit)          unit  = 0;
    যদি (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    যদি (side & 1) nrowa = args.n;

    যদি (args.ldb < MAX(1,args.m)) info = 11;
    যদি (args.lda < MAX(1,nrowa))  info =  9;
    যদি (args.n < 0)               info =  6;
    যদি (args.m < 0)               info =  5;
    যদি (unit < 0)                 info =  4;
    যদি (trans < 0)                info =  3;
    যদি (uplo  < 0)                info =  2;
    যদি (side  < 0)                info =  1;
  }

  যদি (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    যদি (Side == Cblas_বাম)         side  = 1;
    যদি (Side == Cblas_ডান)        side  = 0;

    যদি (Uplo == CblasUpper)        uplo  = 1;
    যদি (Uplo == CblasLower)        uplo  = 0;

    যদি (Trans == CblasNoTrans)     trans = 0;
    যদি (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    যদি (Trans == CblasConjNoTrans) trans = 0;
    যদি (Trans == CblasConjTrans)   trans = 1;
#অন্যথা
    যদি (Trans == CblasConjNoTrans) trans = 2;
    যদি (Trans == CblasConjTrans)   trans = 3;
#endif

    যদি (Diag == CblasUnit)         unit  = 0;
    যদি (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    যদি (side & 1) nrowa = args.n;

    যদি (args.ldb < MAX(1,args.m)) info = 11;
    যদি (args.lda < MAX(1,nrowa))  info =  9;
    যদি (args.n < 0)               info =  6;
    যদি (args.m < 0)               info =  5;
    যদি (unit < 0)                 info =  4;
    যদি (trans < 0)                info =  3;
    যদি (uplo  < 0)                info =  2;
    যদি (side  < 0)                info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

  যদি (args.m == 0 || args.n == 0) প্রত্যাগম;

#যদি defined(TRMM)
#যদি !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#যদি defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
যদি (strcmp(gotoblas_corename(), "armv9sme") == 0
#যদি defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  যদি (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_বাম && m == lda && n == ldb) {
    যদি (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      প্রত্যাগম;
    } অন্যথা যদি (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      প্রত্যাগম;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  যদি ((args.m == 0) || (args.n == 0)) প্রত্যাগম;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
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
  যদি ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  অন্যথা
	args.nthreads = অঙ্ক_cpu_avail(3);
		

  যদি (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, রিক্ত, রিক্ত, sa, sb, 0);

#ifdef SMP
  } অন্যথা {
    যদি (!side) {
      gemm_thread_n(mode, &args, রিক্ত, রিক্ত, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } অন্যথা {
      gemm_thread_m(mode, &args, রিক্ত, রিক্ত, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  প্রত্যাগম;
}

