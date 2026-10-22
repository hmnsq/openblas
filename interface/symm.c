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

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#இல்லையெனில்
#define ERROR_NAME "SSYMM "
#endif
#இல்லையெனில்
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#இல்லையெனில்
#define ERROR_NAME "CSYMM "
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#இல்லையெனில்
#define ERROR_NAME "CHEMM "
#endif
#endif
#இல்லையெனில்
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#இல்லையெனில்
#define ERROR_NAME "CSYMM3M "
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#இல்லையெனில்
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
வெளி எழுத்து* gotoblas_corename(வெறுமை);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#அல்லதெனில் defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#இல்லையெனில்
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#அல்லதெனில் defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#இல்லையெனில்
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

நிலையான முழு (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#இல்லையெனில்
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#இல்லையெனில்
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#இல்லையெனில்
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *SIDE, எழுத்து *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  எழுத்து side_arg  = *SIDE;
  எழுத்து uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#எனில் defined(SMP) && !defined(NO_AFFINITY)
  முழு nodes;
#endif
# எனில் defined(SMP)
  இரட்டை MN;
#endif
  blasint info;
  முழு side;
  முழு uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  எனில் (side_arg  == 'L') side  = 0;
  எனில் (side_arg  == 'R') side  = 1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (வெறுமை *)c;
  args.ldc = *ldC;

  info = 0;

  எனில் (args.ldc < MAX(1, args.m)) info = 12;

  எனில் (!side) {
    args.a = (வெறுமை *)a;
    args.b = (வெறுமை *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    எனில் (args.ldb < MAX(1, args.m)) info =  9;
    எனில் (args.lda < MAX(1, args.m)) info =  7;

  } இல்லையெனில் {
    args.a = (வெறுமை *)b;
    args.b = (வெறுமை *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  எனில் (args.lda < MAX(1, args.m)) info =  9;
  எனில் (args.ldb < MAX(1, args.n)) info =  7;
  }

  எனில் (args.n   < 0)              info =  4;
  எனில் (args.m   < 0)              info =  3;
  எனில் (uplo     < 0)              info =  2;
  எனில் (side     < 0)              info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_SIDE Side, எண்ணிக்கை CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#இல்லையெனில்
	   வெறுமை *valpha,
	   வெறுமை *va, blasint lda,
	   வெறுமை *vb, blasint ldb,
	   வெறுமை *vbeta,
	   வெறுமை *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  முழு side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#எனில் defined(SMP) && !defined(NO_AFFINITY)
  முழு nodes;
#endif
#எனில் defined(SMP)
  இரட்டை MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (வெறுமை *)&alpha;
  args.beta  = (வெறுமை *)&beta;
#இல்லையெனில்
  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;
#endif

  args.c = (வெறுமை *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Side == Cblas_இடது)  side = 0;
    எனில் (Side == Cblas_வலது) side = 1;

    எனில் (Uplo == CblasUpper) uplo  = 0;
    எனில் (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    எனில் (args.ldc < MAX(1, args.m)) info = 12;

    எனில் (!side) {
      args.a = (வெறுமை *)a;
      args.b = (வெறுமை *)b;

      args.lda = lda;
      args.ldb = ldb;

      எனில் (args.ldb < MAX(1, args.m)) info =  9;
      எனில் (args.lda < MAX(1, args.m)) info =  7;

    } இல்லையெனில் {
      args.a = (வெறுமை *)b;
      args.b = (வெறுமை *)a;

      args.lda = ldb;
      args.ldb = lda;

      எனில் (args.lda < MAX(1, args.m)) info =  9;
      எனில் (args.ldb < MAX(1, args.n)) info =  7;
    }

    எனில் (args.n   < 0)              info =  4;
    எனில் (args.m   < 0)              info =  3;
    எனில் (uplo     < 0)              info =  2;
    எனில் (side     < 0)              info =  1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (Side == Cblas_இடது)  side = 1;
    எனில் (Side == Cblas_வலது) side = 0;

    எனில் (Uplo == CblasUpper) uplo  = 1;
    எனில் (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    எனில் (args.ldc < MAX(1, args.m)) info = 12;

    எனில் (!side) {
      args.a = (வெறுமை *)a;
      args.b = (வெறுமை *)b;

      args.lda = lda;
      args.ldb = ldb;

      எனில் (args.ldb < MAX(1, args.m)) info =  9;
      எனில் (args.lda < MAX(1, args.m)) info =  7;

    } இல்லையெனில் {
      args.a = (வெறுமை *)b;
      args.b = (வெறுமை *)a;

      args.lda = ldb;
      args.ldb = lda;

      எனில் (args.lda < MAX(1, args.m)) info =  9;
      எனில் (args.ldb < MAX(1, args.n)) info =  7;
    }

    எனில் (args.n   < 0)              info =  4;
    எனில் (args.m   < 0)              info =  3;
    எனில் (uplo     < 0)              info =  2;
    எனில் (side     < 0)              info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

   எனில் (args.m == 0 || args.n == 0) திரும்பு;
#எனில் !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   எனில் (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     எனில் (Side == Cblas_இடது && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); திரும்பு;
     }
     இல்லையெனில் எனில் (Side == Cblas_இடது && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); திரும்பு;
     }
   }
#endif
#endif

#endif

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
  args.common = வெற்று;
  MN = 2.* (இரட்டை) args.m * (இரட்டை)args.m * (இரட்டை) args.n;
  எனில் (MN <= (SMP_THRESHOLD_MIN * (இரட்டை) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } இல்லையெனில் {
	 args.nthreads = இலக்கம்_cpu_avail(3);
  }
  எனில் (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP
  } இல்லையெனில் {

#ifndef NO_AFFINITY
    nodes = get_இலக்கம்_nodes();

    எனில் (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, வெற்று, வெற்று,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } இல்லையெனில் {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, வெற்று, வெற்று, sa, sb, 0);

#இல்லையெனில்

      GEMM_THREAD(MODE, &args, வெற்று, வெற்று, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  திரும்பு;
}
