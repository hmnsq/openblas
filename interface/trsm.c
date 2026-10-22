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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <ctype.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#ଅନ୍ୟଥା
#define ERROR_NAME "STRSM "
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#ଅନ୍ୟଥା
#define ERROR_NAME "CTRSM "
#endif
#endif
#ଅନ୍ୟଥା
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "STRMM "
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#ଅନ୍ୟଥା
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
ବାହ୍ୟ ଅକ୍ଷର* gotoblas_corename(ଶୂନ୍ୟ);                                
#endif                                                               

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#ଅନ୍ୟଥା
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

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *SIDE, ଅକ୍ଷର *UPLO, ଅକ୍ଷର *TRANS, ଅକ୍ଷର *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  ଅକ୍ଷର side_arg  = *SIDE;
  ଅକ୍ଷର uplo_arg  = *UPLO;
  ଅକ୍ଷର trans_arg = *TRANS;
  ଅକ୍ଷର diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା side;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା unit;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (ଶୂନ୍ୟ *)a;
  args.b = (ଶୂନ୍ୟ *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (ଶୂନ୍ୟ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  ଯଦି (side_arg  == 'L') side  = 0;
  ଯଦି (side_arg  == 'R') side  = 1;

  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
  ଯଦି (trans_arg == 'R') trans = 2;
  ଯଦି (trans_arg == 'C') trans = 3;

  ଯଦି (diag_arg  == 'U') unit  = 0;
  ଯଦି (diag_arg  == 'N') unit  = 1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  ଯଦି (side & 1) nrowa = args.n;

  info = 0;

  ଯଦି (args.ldb < MAX(1,args.m)) info = 11;
  ଯଦି (args.lda < MAX(1,nrowa))  info =  9;
  ଯଦି (args.n < 0)               info =  6;
  ଯଦି (args.m < 0)               info =  5;
  ଯଦି (unit < 0)                 info =  4;
  ଯଦି (trans < 0)                info =  3;
  ଯଦି (uplo  < 0)                info =  2;
  ଯଦି (side  < 0)                info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME)-1);
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order,
	   ଗଣନା CBLAS_SIDE Side,  ଗଣନା CBLAS_UPLO Uplo,
	   ଗଣନା CBLAS_TRANSPOSE Trans, ଗଣନା CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *valpha,
	   ଶୂନ୍ୟ *va, blasint lda,
	   ଶୂନ୍ୟ *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#endif

  blas_arg_t args;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (ଶୂନ୍ୟ *)a;
  args.b = (ଶୂନ୍ୟ *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (ଶୂନ୍ୟ *)&alpha;
#ଅନ୍ୟଥା
  args.beta = (ଶୂନ୍ୟ *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  ଯଦି (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    ଯଦି (Side == Cblas_ବାମ)         side  = 0;
    ଯଦି (Side == Cblas_ଡାହାଣ)        side  = 1;

    ଯଦି (Uplo == CblasUpper)        uplo  = 0;
    ଯଦି (Uplo == CblasLower)        uplo  = 1;

    ଯଦି (Trans == CblasNoTrans)     trans = 0;
    ଯଦି (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    ଯଦି (Trans == CblasConjNoTrans) trans = 0;
    ଯଦି (Trans == CblasConjTrans)   trans = 1;
#ଅନ୍ୟଥା
    ଯଦି (Trans == CblasConjNoTrans) trans = 2;
    ଯଦି (Trans == CblasConjTrans)   trans = 3;
#endif

    ଯଦି (Diag == CblasUnit)          unit  = 0;
    ଯଦି (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    ଯଦି (side & 1) nrowa = args.n;

    ଯଦି (args.ldb < MAX(1,args.m)) info = 11;
    ଯଦି (args.lda < MAX(1,nrowa))  info =  9;
    ଯଦି (args.n < 0)               info =  6;
    ଯଦି (args.m < 0)               info =  5;
    ଯଦି (unit < 0)                 info =  4;
    ଯଦି (trans < 0)                info =  3;
    ଯଦି (uplo  < 0)                info =  2;
    ଯଦି (side  < 0)                info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    ଯଦି (Side == Cblas_ବାମ)         side  = 1;
    ଯଦି (Side == Cblas_ଡାହାଣ)        side  = 0;

    ଯଦି (Uplo == CblasUpper)        uplo  = 1;
    ଯଦି (Uplo == CblasLower)        uplo  = 0;

    ଯଦି (Trans == CblasNoTrans)     trans = 0;
    ଯଦି (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    ଯଦି (Trans == CblasConjNoTrans) trans = 0;
    ଯଦି (Trans == CblasConjTrans)   trans = 1;
#ଅନ୍ୟଥା
    ଯଦି (Trans == CblasConjNoTrans) trans = 2;
    ଯଦି (Trans == CblasConjTrans)   trans = 3;
#endif

    ଯଦି (Diag == CblasUnit)         unit  = 0;
    ଯଦି (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    ଯଦି (side & 1) nrowa = args.n;

    ଯଦି (args.ldb < MAX(1,args.m)) info = 11;
    ଯଦି (args.lda < MAX(1,nrowa))  info =  9;
    ଯଦି (args.n < 0)               info =  6;
    ଯଦି (args.m < 0)               info =  5;
    ଯଦି (unit < 0)                 info =  4;
    ଯଦି (trans < 0)                info =  3;
    ଯଦି (uplo  < 0)                info =  2;
    ଯଦି (side  < 0)                info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

  ଯଦି (args.m == 0 || args.n == 0) ଫେରନ୍ତୁ;

#ଯଦି defined(TRMM)
#ଯଦି !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ଯଦି defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
ଯଦି (strcmp(gotoblas_corename(), "armv9sme") == 0
#ଯଦି defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  ଯଦି (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_ବାମ && m == lda && n == ldb) {
    ଯଦି (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      ଫେରନ୍ତୁ;
    } ଅନ୍ୟଥା ଯଦି (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      ଫେରନ୍ତୁ;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  ଯଦି ((args.m == 0) || (args.n == 0)) ଫେରନ୍ତୁ;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
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
  ଯଦି ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  ଅନ୍ୟଥା
	args.nthreads = ଅଙ୍କ_cpu_avail(3);
		

  ଯଦି (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP
  } ଅନ୍ୟଥା {
    ଯଦି (!side) {
      gemm_thread_n(mode, &args, ଖାଲି, ଖାଲି, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } ଅନ୍ୟଥା {
      gemm_thread_m(mode, &args, ଖାଲି, ଖାଲି, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}

