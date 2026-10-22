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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <ctype.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#இல்லையெனில்
#define ERROR_NAME "STRSM "
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#இல்லையெனில்
#define ERROR_NAME "CTRSM "
#endif
#endif
#இல்லையெனில்
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#இல்லையெனில்
#define ERROR_NAME "STRMM "
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#இல்லையெனில்
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#இல்லையெனில்
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
வெளி எழுத்து* gotoblas_corename(வெறுமை);                                
#endif                                                               

நிலையான முழு (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#இல்லையெனில்
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

வெறுமை NAME(எழுத்து *SIDE, எழுத்து *UPLO, எழுத்து *TRANS, எழுத்து *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  எழுத்து side_arg  = *SIDE;
  எழுத்து uplo_arg  = *UPLO;
  எழுத்து trans_arg = *TRANS;
  எழுத்து diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  முழு side;
  முழு uplo;
  முழு unit;
  முழு trans;
  முழு nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (வெறுமை *)a;
  args.b = (வெறுமை *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (வெறுமை *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  எனில் (side_arg  == 'L') side  = 0;
  எனில் (side_arg  == 'R') side  = 1;

  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
  எனில் (trans_arg == 'R') trans = 2;
  எனில் (trans_arg == 'C') trans = 3;

  எனில் (diag_arg  == 'U') unit  = 0;
  எனில் (diag_arg  == 'N') unit  = 1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  எனில் (side & 1) nrowa = args.n;

  info = 0;

  எனில் (args.ldb < MAX(1,args.m)) info = 11;
  எனில் (args.lda < MAX(1,nrowa))  info =  9;
  எனில் (args.n < 0)               info =  6;
  எனில் (args.m < 0)               info =  5;
  எனில் (unit < 0)                 info =  4;
  எனில் (trans < 0)                info =  3;
  எனில் (uplo  < 0)                info =  2;
  எனில் (side  < 0)                info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME)-1);
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   எண்ணிக்கை CBLAS_SIDE Side,  எண்ணிக்கை CBLAS_UPLO Uplo,
	   எண்ணிக்கை CBLAS_TRANSPOSE Trans, எண்ணிக்கை CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#இல்லையெனில்
	   வெறுமை *valpha,
	   வெறுமை *va, blasint lda,
	   வெறுமை *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#endif

  blas_arg_t args;
  முழு side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (வெறுமை *)a;
  args.b = (வெறுமை *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (வெறுமை *)&alpha;
#இல்லையெனில்
  args.beta = (வெறுமை *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  எனில் (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    எனில் (Side == Cblas_இடது)         side  = 0;
    எனில் (Side == Cblas_வலது)        side  = 1;

    எனில் (Uplo == CblasUpper)        uplo  = 0;
    எனில் (Uplo == CblasLower)        uplo  = 1;

    எனில் (Trans == CblasNoTrans)     trans = 0;
    எனில் (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    எனில் (Trans == CblasConjNoTrans) trans = 0;
    எனில் (Trans == CblasConjTrans)   trans = 1;
#இல்லையெனில்
    எனில் (Trans == CblasConjNoTrans) trans = 2;
    எனில் (Trans == CblasConjTrans)   trans = 3;
#endif

    எனில் (Diag == CblasUnit)          unit  = 0;
    எனில் (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    எனில் (side & 1) nrowa = args.n;

    எனில் (args.ldb < MAX(1,args.m)) info = 11;
    எனில் (args.lda < MAX(1,nrowa))  info =  9;
    எனில் (args.n < 0)               info =  6;
    எனில் (args.m < 0)               info =  5;
    எனில் (unit < 0)                 info =  4;
    எனில் (trans < 0)                info =  3;
    எனில் (uplo  < 0)                info =  2;
    எனில் (side  < 0)                info =  1;
  }

  எனில் (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    எனில் (Side == Cblas_இடது)         side  = 1;
    எனில் (Side == Cblas_வலது)        side  = 0;

    எனில் (Uplo == CblasUpper)        uplo  = 1;
    எனில் (Uplo == CblasLower)        uplo  = 0;

    எனில் (Trans == CblasNoTrans)     trans = 0;
    எனில் (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    எனில் (Trans == CblasConjNoTrans) trans = 0;
    எனில் (Trans == CblasConjTrans)   trans = 1;
#இல்லையெனில்
    எனில் (Trans == CblasConjNoTrans) trans = 2;
    எனில் (Trans == CblasConjTrans)   trans = 3;
#endif

    எனில் (Diag == CblasUnit)         unit  = 0;
    எனில் (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    எனில் (side & 1) nrowa = args.n;

    எனில் (args.ldb < MAX(1,args.m)) info = 11;
    எனில் (args.lda < MAX(1,nrowa))  info =  9;
    எனில் (args.n < 0)               info =  6;
    எனில் (args.m < 0)               info =  5;
    எனில் (unit < 0)                 info =  4;
    எனில் (trans < 0)                info =  3;
    எனில் (uplo  < 0)                info =  2;
    எனில் (side  < 0)                info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

  எனில் (args.m == 0 || args.n == 0) திரும்பு;

#எனில் defined(TRMM)
#எனில் !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  எனில் (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_இடது && m == lda && n == ldb) {
    எனில் (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      திரும்பு;
    } இல்லையெனில் எனில் (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      திரும்பு;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  எனில் ((args.m == 0) || (args.n == 0)) திரும்பு;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
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
  எனில் ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  இல்லையெனில்
	args.nthreads = இலக்கம்_cpu_avail(3);
		

  எனில் (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP
  } இல்லையெனில் {
    எனில் (!side) {
      gemm_thread_n(mode, &args, வெற்று, வெற்று, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } இல்லையெனில் {
      gemm_thread_m(mode, &args, வெற்று, வெற்று, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  திரும்பு;
}

