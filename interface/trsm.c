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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <ctype.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "STRSM "
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CTRSM "
#endif
#endif
#അല്ലെങ്കിൽ
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "STRMM "
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#അല്ലെങ്കിൽ
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
ബാഹ്യ അക്ഷരം* gotoblas_corename(ശൂന്യം);                                
#endif                                                               

സ്ഥിതം പൂർണ്ണസംഖ്യ (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#അല്ലെങ്കിൽ
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

ശൂന്യം NAME(അക്ഷരം *SIDE, അക്ഷരം *UPLO, അക്ഷരം *TRANS, അക്ഷരം *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  അക്ഷരം side_arg  = *SIDE;
  അക്ഷരം uplo_arg  = *UPLO;
  അക്ഷരം trans_arg = *TRANS;
  അക്ഷരം diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  പൂർണ്ണസംഖ്യ side;
  പൂർണ്ണസംഖ്യ uplo;
  പൂർണ്ണസംഖ്യ unit;
  പൂർണ്ണസംഖ്യ trans;
  പൂർണ്ണസംഖ്യ nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (ശൂന്യം *)a;
  args.b = (ശൂന്യം *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (ശൂന്യം *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  എങ്കിൽ (side_arg  == 'L') side  = 0;
  എങ്കിൽ (side_arg  == 'R') side  = 1;

  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'T') trans = 1;
  എങ്കിൽ (trans_arg == 'R') trans = 2;
  എങ്കിൽ (trans_arg == 'C') trans = 3;

  എങ്കിൽ (diag_arg  == 'U') unit  = 0;
  എങ്കിൽ (diag_arg  == 'N') unit  = 1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  എങ്കിൽ (side & 1) nrowa = args.n;

  info = 0;

  എങ്കിൽ (args.ldb < MAX(1,args.m)) info = 11;
  എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  9;
  എങ്കിൽ (args.n < 0)               info =  6;
  എങ്കിൽ (args.m < 0)               info =  5;
  എങ്കിൽ (unit < 0)                 info =  4;
  എങ്കിൽ (trans < 0)                info =  3;
  എങ്കിൽ (uplo  < 0)                info =  2;
  എങ്കിൽ (side  < 0)                info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME)-1);
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order,
	   എണ്ണൽ CBLAS_SIDE Side,  എണ്ണൽ CBLAS_UPLO Uplo,
	   എണ്ണൽ CBLAS_TRANSPOSE Trans, എണ്ണൽ CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#അല്ലെങ്കിൽ
	   ശൂന്യം *valpha,
	   ശൂന്യം *va, blasint lda,
	   ശൂന്യം *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#endif

  blas_arg_t args;
  പൂർണ്ണസംഖ്യ side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (ശൂന്യം *)a;
  args.b = (ശൂന്യം *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (ശൂന്യം *)&alpha;
#അല്ലെങ്കിൽ
  args.beta = (ശൂന്യം *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  എങ്കിൽ (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    എങ്കിൽ (Side == Cblas_ഇടത്)         side  = 0;
    എങ്കിൽ (Side == Cblas_വലത്)        side  = 1;

    എങ്കിൽ (Uplo == CblasUpper)        uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower)        uplo  = 1;

    എങ്കിൽ (Trans == CblasNoTrans)     trans = 0;
    എങ്കിൽ (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 2;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 3;
#endif

    എങ്കിൽ (Diag == CblasUnit)          unit  = 0;
    എങ്കിൽ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    എങ്കിൽ (side & 1) nrowa = args.n;

    എങ്കിൽ (args.ldb < MAX(1,args.m)) info = 11;
    എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  9;
    എങ്കിൽ (args.n < 0)               info =  6;
    എങ്കിൽ (args.m < 0)               info =  5;
    എങ്കിൽ (unit < 0)                 info =  4;
    എങ്കിൽ (trans < 0)                info =  3;
    എങ്കിൽ (uplo  < 0)                info =  2;
    എങ്കിൽ (side  < 0)                info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    എങ്കിൽ (Side == Cblas_ഇടത്)         side  = 1;
    എങ്കിൽ (Side == Cblas_വലത്)        side  = 0;

    എങ്കിൽ (Uplo == CblasUpper)        uplo  = 1;
    എങ്കിൽ (Uplo == CblasLower)        uplo  = 0;

    എങ്കിൽ (Trans == CblasNoTrans)     trans = 0;
    എങ്കിൽ (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 2;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 3;
#endif

    എങ്കിൽ (Diag == CblasUnit)         unit  = 0;
    എങ്കിൽ (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    എങ്കിൽ (side & 1) nrowa = args.n;

    എങ്കിൽ (args.ldb < MAX(1,args.m)) info = 11;
    എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  9;
    എങ്കിൽ (args.n < 0)               info =  6;
    എങ്കിൽ (args.m < 0)               info =  5;
    എങ്കിൽ (unit < 0)                 info =  4;
    എങ്കിൽ (trans < 0)                info =  3;
    എങ്കിൽ (uplo  < 0)                info =  2;
    എങ്കിൽ (side  < 0)                info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

  എങ്കിൽ (args.m == 0 || args.n == 0) തിരികെ;

#എങ്കിൽ defined(TRMM)
#എങ്കിൽ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#എങ്കിൽ defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
എങ്കിൽ (strcmp(gotoblas_corename(), "armv9sme") == 0
#എങ്കിൽ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  എങ്കിൽ (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_ഇടത് && m == lda && n == ldb) {
    എങ്കിൽ (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      തിരികെ;
    } അല്ലെങ്കിൽ എങ്കിൽ (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      തിരികെ;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  എങ്കിൽ ((args.m == 0) || (args.n == 0)) തിരികെ;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
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
  എങ്കിൽ ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  അല്ലെങ്കിൽ
	args.nthreads = അക്കം_cpu_avail(3);
		

  എങ്കിൽ (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP
  } അല്ലെങ്കിൽ {
    എങ്കിൽ (!side) {
      gemm_thread_n(mode, &args, ശൂന്യ, ശൂന്യ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } അല്ലെങ്കിൽ {
      gemm_thread_m(mode, &args, ശൂന്യ, ശൂന്യ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  തിരികെ;
}

