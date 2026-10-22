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

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#ܐܚܪܢܐ
#define ERROR_NAME "SSYMM "
#endif
#ܐܚܪܢܐ
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#ܐܚܪܢܐ
#define ERROR_NAME "CSYMM "
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#ܐܚܪܢܐ
#define ERROR_NAME "CHEMM "
#endif
#endif
#ܐܚܪܢܐ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#ܐܚܪܢܐ
#define ERROR_NAME "CSYMM3M "
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#ܐܚܪܢܐ
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#وإلاإذا defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#ܐܚܪܢܐ
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ܐܚܪܢܐ
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#ܐܚܪܢܐ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#ܐܚܪܢܐ
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#ܐܚܪܢܐ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ܐܢ defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#endif
# ܐܢ defined(SMP)
  مزدوج MN;
#endif
  blasint info;
  صحيح side;
  صحيح uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ܐܢ (side_arg  == 'L') side  = 0;
  ܐܢ (side_arg  == 'R') side  = 1;

  ܐܢ (uplo_arg  == 'U') uplo  = 0;
  ܐܢ (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (فراغ *)c;
  args.ldc = *ldC;

  info = 0;

  ܐܢ (args.ldc < MAX(1, args.m)) info = 12;

  ܐܢ (!side) {
    args.a = (فراغ *)a;
    args.b = (فراغ *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    ܐܢ (args.ldb < MAX(1, args.m)) info =  9;
    ܐܢ (args.lda < MAX(1, args.m)) info =  7;

  } ܐܚܪܢܐ {
    args.a = (فراغ *)b;
    args.b = (فراغ *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  ܐܢ (args.lda < MAX(1, args.m)) info =  9;
  ܐܢ (args.ldb < MAX(1, args.n)) info =  7;
  }

  ܐܢ (args.n   < 0)              info =  4;
  ܐܢ (args.m   < 0)              info =  3;
  ܐܢ (uplo     < 0)              info =  2;
  ܐܢ (side     < 0)              info =  1;

  ܐܢ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ܐܚܪܢܐ
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
	   فراغ *vbeta,
	   فراغ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  صحيح side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ܐܢ defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#endif
#ܐܢ defined(SMP)
  مزدوج MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (فراغ *)&alpha;
  args.beta  = (فراغ *)&beta;
#ܐܚܪܢܐ
  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;
#endif

  args.c = (فراغ *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ܐܢ (order == CblasColMajor) {
    ܐܢ (Side == CblasLeft)  side = 0;
    ܐܢ (Side == CblasRight) side = 1;

    ܐܢ (Uplo == CblasUpper) uplo  = 0;
    ܐܢ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    ܐܢ (args.ldc < MAX(1, args.m)) info = 12;

    ܐܢ (!side) {
      args.a = (فراغ *)a;
      args.b = (فراغ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ܐܢ (args.ldb < MAX(1, args.m)) info =  9;
      ܐܢ (args.lda < MAX(1, args.m)) info =  7;

    } ܐܚܪܢܐ {
      args.a = (فراغ *)b;
      args.b = (فراغ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ܐܢ (args.lda < MAX(1, args.m)) info =  9;
      ܐܢ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ܐܢ (args.n   < 0)              info =  4;
    ܐܢ (args.m   < 0)              info =  3;
    ܐܢ (uplo     < 0)              info =  2;
    ܐܢ (side     < 0)              info =  1;
  }

  ܐܢ (order == CblasRowMajor) {
    ܐܢ (Side == CblasLeft)  side = 1;
    ܐܢ (Side == CblasRight) side = 0;

    ܐܢ (Uplo == CblasUpper) uplo  = 1;
    ܐܢ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    ܐܢ (args.ldc < MAX(1, args.m)) info = 12;

    ܐܢ (!side) {
      args.a = (فراغ *)a;
      args.b = (فراغ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ܐܢ (args.ldb < MAX(1, args.m)) info =  9;
      ܐܢ (args.lda < MAX(1, args.m)) info =  7;

    } ܐܚܪܢܐ {
      args.a = (فراغ *)b;
      args.b = (فراغ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ܐܢ (args.lda < MAX(1, args.m)) info =  9;
      ܐܢ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ܐܢ (args.n   < 0)              info =  4;
    ܐܢ (args.m   < 0)              info =  3;
    ܐܢ (uplo     < 0)              info =  2;
    ܐܢ (side     < 0)              info =  1;
  }

  ܐܢ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

   ܐܢ (args.m == 0 || args.n == 0) ܦܢܐ;
#ܐܢ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ܐܢ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ܐܢ defined(DYNAMIC_ARCH)
ܐܢ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ܐܢ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ܐܢ (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     ܐܢ (Side == CblasLeft && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ܦܢܐ;
     }
     ܐܚܪܢܐ ܐܢ (Side == CblasLeft && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ܦܢܐ;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ܐܢ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = عدم;
  MN = 2.* (مزدوج) args.m * (مزدوج)args.m * (مزدوج) args.n;
  ܐܢ (MN <= (SMP_THRESHOLD_MIN * (مزدوج) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } ܐܚܪܢܐ {
	 args.nthreads = num_cpu_avail(3);
  }
  ܐܢ (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } ܐܚܪܢܐ {

#ifndef NO_AFFINITY
    nodes = get_num_nodes();

    ܐܢ (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ܐܚܪܢܐ {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, عدم, عدم, sa, sb, 0);

#ܐܚܪܢܐ

      GEMM_THREAD(MODE, &args, عدم, عدم, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  ܦܢܐ;
}
