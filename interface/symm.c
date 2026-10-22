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

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "SSYMM "
#endif
#ଅନ୍ୟଥା
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "CSYMM "
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#ଅନ୍ୟଥା
#define ERROR_NAME "CHEMM "
#endif
#endif
#ଅନ୍ୟଥା
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#ଅନ୍ୟଥା
#define ERROR_NAME "CSYMM3M "
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#ଅନ୍ୟଥା
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
ବାହ୍ୟ ଅକ୍ଷର* gotoblas_corename(ଶୂନ୍ୟ);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#ନଚେତଯଦି defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#ଅନ୍ୟଥା
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#ନଚେତଯଦି defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ଅନ୍ୟଥା
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#ଅନ୍ୟଥା
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#ଅନ୍ୟଥା
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#ଅନ୍ୟଥା
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ଯଦି defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *SIDE, ଅକ୍ଷର *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ଅକ୍ଷର side_arg  = *SIDE;
  ଅକ୍ଷର uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ଯଦି defined(SMP) && !defined(NO_AFFINITY)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nodes;
#endif
# ଯଦି defined(SMP)
  ଦ୍ୱିଗୁଣ MN;
#endif
  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା side;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (ଶୂନ୍ୟ *)alpha;
  args.beta  = (ଶୂନ୍ୟ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ଯଦି (side_arg  == 'L') side  = 0;
  ଯଦି (side_arg  == 'R') side  = 1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (ଶୂନ୍ୟ *)c;
  args.ldc = *ldC;

  info = 0;

  ଯଦି (args.ldc < MAX(1, args.m)) info = 12;

  ଯଦି (!side) {
    args.a = (ଶୂନ୍ୟ *)a;
    args.b = (ଶୂନ୍ୟ *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    ଯଦି (args.ldb < MAX(1, args.m)) info =  9;
    ଯଦି (args.lda < MAX(1, args.m)) info =  7;

  } ଅନ୍ୟଥା {
    args.a = (ଶୂନ୍ୟ *)b;
    args.b = (ଶୂନ୍ୟ *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  ଯଦି (args.lda < MAX(1, args.m)) info =  9;
  ଯଦି (args.ldb < MAX(1, args.n)) info =  7;
  }

  ଯଦି (args.n   < 0)              info =  4;
  ଯଦି (args.m   < 0)              info =  3;
  ଯଦି (uplo     < 0)              info =  2;
  ଯଦି (side     < 0)              info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_SIDE Side, ଗଣନା CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ *valpha,
	   ଶୂନ୍ୟ *va, blasint lda,
	   ଶୂନ୍ୟ *vb, blasint ldb,
	   ଶୂନ୍ୟ *vbeta,
	   ଶୂନ୍ୟ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ଯଦି defined(SMP) && !defined(NO_AFFINITY)
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nodes;
#endif
#ଯଦି defined(SMP)
  ଦ୍ୱିଗୁଣ MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ଶୂନ୍ୟ *)&alpha;
  args.beta  = (ଶୂନ୍ୟ *)&beta;
#ଅନ୍ୟଥା
  args.alpha = (ଶୂନ୍ୟ *)alpha;
  args.beta  = (ଶୂନ୍ୟ *)beta;
#endif

  args.c = (ଶୂନ୍ୟ *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (Side == Cblas_ବାମ)  side = 0;
    ଯଦି (Side == Cblas_ଡାହାଣ) side = 1;

    ଯଦି (Uplo == CblasUpper) uplo  = 0;
    ଯଦି (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    ଯଦି (args.ldc < MAX(1, args.m)) info = 12;

    ଯଦି (!side) {
      args.a = (ଶୂନ୍ୟ *)a;
      args.b = (ଶୂନ୍ୟ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ଯଦି (args.ldb < MAX(1, args.m)) info =  9;
      ଯଦି (args.lda < MAX(1, args.m)) info =  7;

    } ଅନ୍ୟଥା {
      args.a = (ଶୂନ୍ୟ *)b;
      args.b = (ଶୂନ୍ୟ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ଯଦି (args.lda < MAX(1, args.m)) info =  9;
      ଯଦି (args.ldb < MAX(1, args.n)) info =  7;
    }

    ଯଦି (args.n   < 0)              info =  4;
    ଯଦି (args.m   < 0)              info =  3;
    ଯଦି (uplo     < 0)              info =  2;
    ଯଦି (side     < 0)              info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (Side == Cblas_ବାମ)  side = 1;
    ଯଦି (Side == Cblas_ଡାହାଣ) side = 0;

    ଯଦି (Uplo == CblasUpper) uplo  = 1;
    ଯଦି (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    ଯଦି (args.ldc < MAX(1, args.m)) info = 12;

    ଯଦି (!side) {
      args.a = (ଶୂନ୍ୟ *)a;
      args.b = (ଶୂନ୍ୟ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ଯଦି (args.ldb < MAX(1, args.m)) info =  9;
      ଯଦି (args.lda < MAX(1, args.m)) info =  7;

    } ଅନ୍ୟଥା {
      args.a = (ଶୂନ୍ୟ *)b;
      args.b = (ଶୂନ୍ୟ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ଯଦି (args.lda < MAX(1, args.m)) info =  9;
      ଯଦି (args.ldb < MAX(1, args.n)) info =  7;
    }

    ଯଦି (args.n   < 0)              info =  4;
    ଯଦି (args.m   < 0)              info =  3;
    ଯଦି (uplo     < 0)              info =  2;
    ଯଦି (side     < 0)              info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

   ଯଦି (args.m == 0 || args.n == 0) ଫେରନ୍ତୁ;
#ଯଦି !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ଯଦି defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ଯଦି defined(DYNAMIC_ARCH)
ଯଦି (strcmp(gotoblas_corename(), "armv9sme") == 0
#ଯଦି defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ଯଦି (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     ଯଦି (Side == Cblas_ବାମ && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ଫେରନ୍ତୁ;
     }
     ଅନ୍ୟଥା ଯଦି (Side == Cblas_ବାମ && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ଫେରନ୍ତୁ;
     }
   }
#endif
#endif

#endif

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
  args.common = ଖାଲି;
  MN = 2.* (ଦ୍ୱିଗୁଣ) args.m * (ଦ୍ୱିଗୁଣ)args.m * (ଦ୍ୱିଗୁଣ) args.n;
  ଯଦି (MN <= (SMP_THRESHOLD_MIN * (ଦ୍ୱିଗୁଣ) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } ଅନ୍ୟଥା {
	 args.nthreads = ଅଙ୍କ_cpu_avail(3);
  }
  ଯଦି (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP
  } ଅନ୍ୟଥା {

#ifndef NO_AFFINITY
    nodes = get_ଅଙ୍କ_nodes();

    ଯଦି (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, ଖାଲି, ଖାଲି,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ଅନ୍ୟଥା {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ଅନ୍ୟଥା

      GEMM_THREAD(MODE, &args, ଖାଲି, ଖାଲି, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  ଫେରନ୍ତୁ;
}
