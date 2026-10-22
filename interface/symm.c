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
#وإلا
#define ERROR_NAME "SSYMM "
#endif
#وإلا
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#وإلا
#define ERROR_NAME "CSYMM "
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#وإلا
#define ERROR_NAME "CHEMM "
#endif
#endif
#وإلا
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#وإلا
#define ERROR_NAME "CSYMM3M "
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#وإلا
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
#وإلا
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#وإلا
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#وإلا
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#وإلا
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#وإلا
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#وإلا
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
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

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#endif
# إذا defined(SMP)
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

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (فراغ *)c;
  args.ldc = *ldC;

  info = 0;

  إذا (args.ldc < MAX(1, args.m)) info = 12;

  إذا (!side) {
    args.a = (فراغ *)a;
    args.b = (فراغ *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    إذا (args.ldb < MAX(1, args.m)) info =  9;
    إذا (args.lda < MAX(1, args.m)) info =  7;

  } وإلا {
    args.a = (فراغ *)b;
    args.b = (فراغ *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  إذا (args.lda < MAX(1, args.m)) info =  9;
  إذا (args.ldb < MAX(1, args.n)) info =  7;
  }

  إذا (args.n   < 0)              info =  4;
  إذا (args.m   < 0)              info =  3;
  إذا (uplo     < 0)              info =  2;
  إذا (side     < 0)              info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#وإلا
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

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#endif
#إذا defined(SMP)
  مزدوج MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (فراغ *)&alpha;
  args.beta  = (فراغ *)&beta;
#وإلا
  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;
#endif

  args.c = (فراغ *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  إذا (order == CblasColMajor) {
    إذا (Side == Cblas_کيري)  side = 0;
    إذا (Side == Cblas_کانن) side = 1;

    إذا (Uplo == CblasUpper) uplo  = 0;
    إذا (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    إذا (args.ldc < MAX(1, args.m)) info = 12;

    إذا (!side) {
      args.a = (فراغ *)a;
      args.b = (فراغ *)b;

      args.lda = lda;
      args.ldb = ldb;

      إذا (args.ldb < MAX(1, args.m)) info =  9;
      إذا (args.lda < MAX(1, args.m)) info =  7;

    } وإلا {
      args.a = (فراغ *)b;
      args.b = (فراغ *)a;

      args.lda = ldb;
      args.ldb = lda;

      إذا (args.lda < MAX(1, args.m)) info =  9;
      إذا (args.ldb < MAX(1, args.n)) info =  7;
    }

    إذا (args.n   < 0)              info =  4;
    إذا (args.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (order == CblasRowMajor) {
    إذا (Side == Cblas_کيري)  side = 1;
    إذا (Side == Cblas_کانن) side = 0;

    إذا (Uplo == CblasUpper) uplo  = 1;
    إذا (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    إذا (args.ldc < MAX(1, args.m)) info = 12;

    إذا (!side) {
      args.a = (فراغ *)a;
      args.b = (فراغ *)b;

      args.lda = lda;
      args.ldb = ldb;

      إذا (args.ldb < MAX(1, args.m)) info =  9;
      إذا (args.lda < MAX(1, args.m)) info =  7;

    } وإلا {
      args.a = (فراغ *)b;
      args.b = (فراغ *)a;

      args.lda = ldb;
      args.ldb = lda;

      إذا (args.lda < MAX(1, args.m)) info =  9;
      إذا (args.ldb < MAX(1, args.n)) info =  7;
    }

    إذا (args.n   < 0)              info =  4;
    إذا (args.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

   إذا (args.m == 0 || args.n == 0) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   إذا (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     إذا (Side == Cblas_کيري && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
     وإلا إذا (Side == Cblas_کيري && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
   }
#endif
#endif

#endif

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
  args.common = عدم;
  MN = 2.* (مزدوج) args.m * (مزدوج)args.m * (مزدوج) args.n;
  إذا (MN <= (SMP_THRESHOLD_MIN * (مزدوج) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } وإلا {
	 args.nthreads = اڠک_cpu_avail(3);
  }
  إذا (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } وإلا {

#ifndef NO_AFFINITY
    nodes = get_اڠک_nodes();

    إذا (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } وإلا {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, عدم, عدم, sa, sb, 0);

#وإلا

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

  أرجع;
}
