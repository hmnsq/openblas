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

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "SSYMM "
#endif
#അല്ലെങ്കിൽ
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CSYMM "
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CHEMM "
#endif
#endif
#അല്ലെങ്കിൽ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CSYMM3M "
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
ബാഹ്യ അക്ഷരം* gotoblas_corename(ശൂന്യം);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#അല്ലെങ്കിൽ
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#അല്ലെങ്കിൽ
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#അല്ലെങ്കിൽ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#അല്ലെങ്കിൽ
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#അല്ലെങ്കിൽ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#എങ്കിൽ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *SIDE, അക്ഷരം *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  അക്ഷരം side_arg  = *SIDE;
  അക്ഷരം uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#എങ്കിൽ defined(SMP) && !defined(NO_AFFINITY)
  പൂർണ്ണസംഖ്യ nodes;
#endif
# എങ്കിൽ defined(SMP)
  ഇരട്ട MN;
#endif
  blasint info;
  പൂർണ്ണസംഖ്യ side;
  പൂർണ്ണസംഖ്യ uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (ശൂന്യം *)alpha;
  args.beta  = (ശൂന്യം *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  എങ്കിൽ (side_arg  == 'L') side  = 0;
  എങ്കിൽ (side_arg  == 'R') side  = 1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (ശൂന്യം *)c;
  args.ldc = *ldC;

  info = 0;

  എങ്കിൽ (args.ldc < MAX(1, args.m)) info = 12;

  എങ്കിൽ (!side) {
    args.a = (ശൂന്യം *)a;
    args.b = (ശൂന്യം *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    എങ്കിൽ (args.ldb < MAX(1, args.m)) info =  9;
    എങ്കിൽ (args.lda < MAX(1, args.m)) info =  7;

  } അല്ലെങ്കിൽ {
    args.a = (ശൂന്യം *)b;
    args.b = (ശൂന്യം *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  എങ്കിൽ (args.lda < MAX(1, args.m)) info =  9;
  എങ്കിൽ (args.ldb < MAX(1, args.n)) info =  7;
  }

  എങ്കിൽ (args.n   < 0)              info =  4;
  എങ്കിൽ (args.m   < 0)              info =  3;
  എങ്കിൽ (uplo     < 0)              info =  2;
  എങ്കിൽ (side     < 0)              info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_SIDE Side, എണ്ണൽ CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#അല്ലെങ്കിൽ
	   ശൂന്യം *valpha,
	   ശൂന്യം *va, blasint lda,
	   ശൂന്യം *vb, blasint ldb,
	   ശൂന്യം *vbeta,
	   ശൂന്യം *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  പൂർണ്ണസംഖ്യ side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#എങ്കിൽ defined(SMP) && !defined(NO_AFFINITY)
  പൂർണ്ണസംഖ്യ nodes;
#endif
#എങ്കിൽ defined(SMP)
  ഇരട്ട MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ശൂന്യം *)&alpha;
  args.beta  = (ശൂന്യം *)&beta;
#അല്ലെങ്കിൽ
  args.alpha = (ശൂന്യം *)alpha;
  args.beta  = (ശൂന്യം *)beta;
#endif

  args.c = (ശൂന്യം *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (Side == Cblas_ഇടത്)  side = 0;
    എങ്കിൽ (Side == Cblas_വലത്) side = 1;

    എങ്കിൽ (Uplo == CblasUpper) uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    എങ്കിൽ (args.ldc < MAX(1, args.m)) info = 12;

    എങ്കിൽ (!side) {
      args.a = (ശൂന്യം *)a;
      args.b = (ശൂന്യം *)b;

      args.lda = lda;
      args.ldb = ldb;

      എങ്കിൽ (args.ldb < MAX(1, args.m)) info =  9;
      എങ്കിൽ (args.lda < MAX(1, args.m)) info =  7;

    } അല്ലെങ്കിൽ {
      args.a = (ശൂന്യം *)b;
      args.b = (ശൂന്യം *)a;

      args.lda = ldb;
      args.ldb = lda;

      എങ്കിൽ (args.lda < MAX(1, args.m)) info =  9;
      എങ്കിൽ (args.ldb < MAX(1, args.n)) info =  7;
    }

    എങ്കിൽ (args.n   < 0)              info =  4;
    എങ്കിൽ (args.m   < 0)              info =  3;
    എങ്കിൽ (uplo     < 0)              info =  2;
    എങ്കിൽ (side     < 0)              info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (Side == Cblas_ഇടത്)  side = 1;
    എങ്കിൽ (Side == Cblas_വലത്) side = 0;

    എങ്കിൽ (Uplo == CblasUpper) uplo  = 1;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    എങ്കിൽ (args.ldc < MAX(1, args.m)) info = 12;

    എങ്കിൽ (!side) {
      args.a = (ശൂന്യം *)a;
      args.b = (ശൂന്യം *)b;

      args.lda = lda;
      args.ldb = ldb;

      എങ്കിൽ (args.ldb < MAX(1, args.m)) info =  9;
      എങ്കിൽ (args.lda < MAX(1, args.m)) info =  7;

    } അല്ലെങ്കിൽ {
      args.a = (ശൂന്യം *)b;
      args.b = (ശൂന്യം *)a;

      args.lda = ldb;
      args.ldb = lda;

      എങ്കിൽ (args.lda < MAX(1, args.m)) info =  9;
      എങ്കിൽ (args.ldb < MAX(1, args.n)) info =  7;
    }

    എങ്കിൽ (args.n   < 0)              info =  4;
    എങ്കിൽ (args.m   < 0)              info =  3;
    എങ്കിൽ (uplo     < 0)              info =  2;
    എങ്കിൽ (side     < 0)              info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

   എങ്കിൽ (args.m == 0 || args.n == 0) തിരികെ;
#എങ്കിൽ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#എങ്കിൽ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
എങ്കിൽ (strcmp(gotoblas_corename(), "armv9sme") == 0
#എങ്കിൽ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   എങ്കിൽ (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     എങ്കിൽ (Side == Cblas_ഇടത് && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); തിരികെ;
     }
     അല്ലെങ്കിൽ എങ്കിൽ (Side == Cblas_ഇടത് && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); തിരികെ;
     }
   }
#endif
#endif

#endif

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
  args.common = ശൂന്യ;
  MN = 2.* (ഇരട്ട) args.m * (ഇരട്ട)args.m * (ഇരട്ട) args.n;
  എങ്കിൽ (MN <= (SMP_THRESHOLD_MIN * (ഇരട്ട) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } അല്ലെങ്കിൽ {
	 args.nthreads = അക്കം_cpu_avail(3);
  }
  എങ്കിൽ (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP
  } അല്ലെങ്കിൽ {

#ifndef NO_AFFINITY
    nodes = get_അക്കം_nodes();

    എങ്കിൽ (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, ശൂന്യ, ശൂന്യ,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } അല്ലെങ്കിൽ {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#അല്ലെങ്കിൽ

      GEMM_THREAD(MODE, &args, ശൂന്യ, ശൂന്യ, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

#endif

#ifndef NO_AFFINITY
    }
#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side)? args.m * (args.m / 2 + args.n) : args.n * (args.m + args.n / 2),
		       (!side)? 2 * args.m * args.m * args.n : 2 * args.m * args.n * args.n);

  IDEBUG_END;

  തിരികെ;
}
