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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SSYMM "
#endif
#ಅನ್ಯಥಾ
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CSYMM "
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHEMM "
#endif
#endif
#ಅನ್ಯಥಾ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CSYMM3M "
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
ಬಾಹ್ಯ ವರ್ಣ* gotoblas_corename(ಶೂನ್ಯ);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#ಅನ್ಯಥಾ
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ಅನ್ಯಥಾ
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#ಅನ್ಯಥಾ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#ಅನ್ಯಥಾ
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#ಅನ್ಯಥಾ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *SIDE, ವರ್ಣ *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ವರ್ಣ side_arg  = *SIDE;
  ವರ್ಣ uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ಯದಿ defined(SMP) && !defined(NO_AFFINITY)
  ಪೂರ್ಣಾಙ್ಕ nodes;
#endif
# ಯದಿ defined(SMP)
  ದ್ವಿಗುಣ MN;
#endif
  blasint info;
  ಪೂರ್ಣಾಙ್ಕ side;
  ಪೂರ್ಣಾಙ್ಕ uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ಯದಿ (side_arg  == 'L') side  = 0;
  ಯದಿ (side_arg  == 'R') side  = 1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (ಶೂನ್ಯ *)c;
  args.ldc = *ldC;

  info = 0;

  ಯದಿ (args.ldc < MAX(1, args.m)) info = 12;

  ಯದಿ (!side) {
    args.a = (ಶೂನ್ಯ *)a;
    args.b = (ಶೂನ್ಯ *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    ಯದಿ (args.ldb < MAX(1, args.m)) info =  9;
    ಯದಿ (args.lda < MAX(1, args.m)) info =  7;

  } ಅನ್ಯಥಾ {
    args.a = (ಶೂನ್ಯ *)b;
    args.b = (ಶೂನ್ಯ *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  ಯದಿ (args.lda < MAX(1, args.m)) info =  9;
  ಯದಿ (args.ldb < MAX(1, args.n)) info =  7;
  }

  ಯದಿ (args.n   < 0)              info =  4;
  ಯದಿ (args.m   < 0)              info =  3;
  ಯದಿ (uplo     < 0)              info =  2;
  ಯದಿ (side     < 0)              info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_SIDE Side, ಗಣನಾ CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb,
	   ಶೂನ್ಯ *vbeta,
	   ಶೂನ್ಯ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  ಪೂರ್ಣಾಙ್ಕ side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ಯದಿ defined(SMP) && !defined(NO_AFFINITY)
  ಪೂರ್ಣಾಙ್ಕ nodes;
#endif
#ಯದಿ defined(SMP)
  ದ್ವಿಗುಣ MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (ಶೂನ್ಯ *)&alpha;
  args.beta  = (ಶೂನ್ಯ *)&beta;
#ಅನ್ಯಥಾ
  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;
#endif

  args.c = (ಶೂನ್ಯ *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Side == Cblas_ಎಡ)  side = 0;
    ಯದಿ (Side == Cblas_ಬಲ) side = 1;

    ಯದಿ (Uplo == CblasUpper) uplo  = 0;
    ಯದಿ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    ಯದಿ (args.ldc < MAX(1, args.m)) info = 12;

    ಯದಿ (!side) {
      args.a = (ಶೂನ್ಯ *)a;
      args.b = (ಶೂನ್ಯ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ಯದಿ (args.ldb < MAX(1, args.m)) info =  9;
      ಯದಿ (args.lda < MAX(1, args.m)) info =  7;

    } ಅನ್ಯಥಾ {
      args.a = (ಶೂನ್ಯ *)b;
      args.b = (ಶೂನ್ಯ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ಯದಿ (args.lda < MAX(1, args.m)) info =  9;
      ಯದಿ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ಯದಿ (args.n   < 0)              info =  4;
    ಯದಿ (args.m   < 0)              info =  3;
    ಯದಿ (uplo     < 0)              info =  2;
    ಯದಿ (side     < 0)              info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Side == Cblas_ಎಡ)  side = 1;
    ಯದಿ (Side == Cblas_ಬಲ) side = 0;

    ಯದಿ (Uplo == CblasUpper) uplo  = 1;
    ಯದಿ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    ಯದಿ (args.ldc < MAX(1, args.m)) info = 12;

    ಯದಿ (!side) {
      args.a = (ಶೂನ್ಯ *)a;
      args.b = (ಶೂನ್ಯ *)b;

      args.lda = lda;
      args.ldb = ldb;

      ಯದಿ (args.ldb < MAX(1, args.m)) info =  9;
      ಯದಿ (args.lda < MAX(1, args.m)) info =  7;

    } ಅನ್ಯಥಾ {
      args.a = (ಶೂನ್ಯ *)b;
      args.b = (ಶೂನ್ಯ *)a;

      args.lda = ldb;
      args.ldb = lda;

      ಯದಿ (args.lda < MAX(1, args.m)) info =  9;
      ಯದಿ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ಯದಿ (args.n   < 0)              info =  4;
    ಯದಿ (args.m   < 0)              info =  3;
    ಯದಿ (uplo     < 0)              info =  2;
    ಯದಿ (side     < 0)              info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

   ಯದಿ (args.m == 0 || args.n == 0) ಪ್ರತ್ಯಾಗಮ;
#ಯದಿ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಯದಿ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
ಯದಿ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಯದಿ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ಯದಿ (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     ಯದಿ (Side == Cblas_ಎಡ && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ಪ್ರತ್ಯಾಗಮ;
     }
     ಅನ್ಯಥಾ ಯದಿ (Side == Cblas_ಎಡ && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ಪ್ರತ್ಯಾಗಮ;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = ರಿಕ್ತ;
  MN = 2.* (ದ್ವಿಗುಣ) args.m * (ದ್ವಿಗುಣ)args.m * (ದ್ವಿಗುಣ) args.n;
  ಯದಿ (MN <= (SMP_THRESHOLD_MIN * (ದ್ವಿಗುಣ) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } ಅನ್ಯಥಾ {
	 args.nthreads = ಅಂಕೆ_cpu_avail(3);
  }
  ಯದಿ (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP
  } ಅನ್ಯಥಾ {

#ifndef NO_AFFINITY
    nodes = get_ಅಂಕೆ_nodes();

    ಯದಿ (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, ರಿಕ್ತ, ರಿಕ್ತ,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ಅನ್ಯಥಾ {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ಅನ್ಯಥಾ

      GEMM_THREAD(MODE, &args, ರಿಕ್ತ, ರಿಕ್ತ, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  ಪ್ರತ್ಯಾಗಮ;
}
