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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "SSYMM "
#endif
#ບໍ່ດັ່ງນັ້ນ
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "CSYMM "
#endif
#ບໍ່ດັ່ງນັ້ນ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "CHEMM "
#endif
#endif
#ບໍ່ດັ່ງນັ້ນ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "CSYMM3M "
#endif
#ບໍ່ດັ່ງນັ້ນ
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
extern var_na* gotoblas_corename(void);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#elif defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#ບໍ່ດັ່ງນັ້ນ
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#ບໍ່ດັ່ງນັ້ນ
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#elif defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ບໍ່ດັ່ງນັ້ນ
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

static ຈຳນວນເຕັມ (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ຖ້າ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#ບໍ່ດັ່ງນັ້ນ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ຖ້າ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#ບໍ່ດັ່ງນັ້ນ
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ຖ້າ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#ບໍ່ດັ່ງນັ້ນ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ຖ້າ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *SIDE, var_na *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  var_na side_arg  = *SIDE;
  var_na uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ຖ້າ defined(SMP) && !defined(NO_AFFINITY)
  ຈຳນວນເຕັມ nodes;
#endif
# ຖ້າ defined(SMP)
  dviguu_na MN;
#endif
  blasint info;
  ຈຳນວນເຕັມ side;
  ຈຳນວນເຕັມ uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (void *)alpha;
  args.beta  = (void *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ຖ້າ (side_arg  == 'L') side  = 0;
  ຖ້າ (side_arg  == 'R') side  = 1;

  ຖ້າ (uplo_arg  == 'U') uplo  = 0;
  ຖ້າ (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (void *)c;
  args.ldc = *ldC;

  info = 0;

  ຖ້າ (args.ldc < MAX(1, args.m)) info = 12;

  ຖ້າ (!side) {
    args.a = (void *)a;
    args.b = (void *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    ຖ້າ (args.ldb < MAX(1, args.m)) info =  9;
    ຖ້າ (args.lda < MAX(1, args.m)) info =  7;

  } ບໍ່ດັ່ງນັ້ນ {
    args.a = (void *)b;
    args.b = (void *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  ຖ້າ (args.lda < MAX(1, args.m)) info =  9;
  ຖ້າ (args.ldb < MAX(1, args.n)) info =  7;
  }

  ຖ້າ (args.n   < 0)              info =  4;
  ຖ້າ (args.m   < 0)              info =  3;
  ຖ້າ (uplo     < 0)              info =  2;
  ຖ້າ (side     < 0)              info =  1;

  ຖ້າ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#ບໍ່ດັ່ງນັ້ນ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_SIDE Side, ga_nanaa CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ບໍ່ດັ່ງນັ້ນ
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb,
	   void *vbeta,
	   void *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  ຈຳນວນເຕັມ side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ຖ້າ defined(SMP) && !defined(NO_AFFINITY)
  ຈຳນວນເຕັມ nodes;
#endif
#ຖ້າ defined(SMP)
  dviguu_na MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (void *)&alpha;
  args.beta  = (void *)&beta;
#ບໍ່ດັ່ງນັ້ນ
  args.alpha = (void *)alpha;
  args.beta  = (void *)beta;
#endif

  args.c = (void *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ຖ້າ (order == CblasColMajor) {
    ຖ້າ (Side == CblasLeft)  side = 0;
    ຖ້າ (Side == CblasRight) side = 1;

    ຖ້າ (Uplo == CblasUpper) uplo  = 0;
    ຖ້າ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    ຖ້າ (args.ldc < MAX(1, args.m)) info = 12;

    ຖ້າ (!side) {
      args.a = (void *)a;
      args.b = (void *)b;

      args.lda = lda;
      args.ldb = ldb;

      ຖ້າ (args.ldb < MAX(1, args.m)) info =  9;
      ຖ້າ (args.lda < MAX(1, args.m)) info =  7;

    } ບໍ່ດັ່ງນັ້ນ {
      args.a = (void *)b;
      args.b = (void *)a;

      args.lda = ldb;
      args.ldb = lda;

      ຖ້າ (args.lda < MAX(1, args.m)) info =  9;
      ຖ້າ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ຖ້າ (args.n   < 0)              info =  4;
    ຖ້າ (args.m   < 0)              info =  3;
    ຖ້າ (uplo     < 0)              info =  2;
    ຖ້າ (side     < 0)              info =  1;
  }

  ຖ້າ (order == CblasRowMajor) {
    ຖ້າ (Side == CblasLeft)  side = 1;
    ຖ້າ (Side == CblasRight) side = 0;

    ຖ້າ (Uplo == CblasUpper) uplo  = 1;
    ຖ້າ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    ຖ້າ (args.ldc < MAX(1, args.m)) info = 12;

    ຖ້າ (!side) {
      args.a = (void *)a;
      args.b = (void *)b;

      args.lda = lda;
      args.ldb = ldb;

      ຖ້າ (args.ldb < MAX(1, args.m)) info =  9;
      ຖ້າ (args.lda < MAX(1, args.m)) info =  7;

    } ບໍ່ດັ່ງນັ້ນ {
      args.a = (void *)b;
      args.b = (void *)a;

      args.lda = ldb;
      args.ldb = lda;

      ຖ້າ (args.lda < MAX(1, args.m)) info =  9;
      ຖ້າ (args.ldb < MAX(1, args.n)) info =  7;
    }

    ຖ້າ (args.n   < 0)              info =  4;
    ຖ້າ (args.m   < 0)              info =  3;
    ຖ້າ (uplo     < 0)              info =  2;
    ຖ້າ (side     < 0)              info =  1;
  }

  ຖ້າ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

   ຖ້າ (args.m == 0 || args.n == 0) ສົ່ງຄືນ;
#ຖ້າ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ຖ້າ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ຖ້າ defined(DYNAMIC_ARCH)
ຖ້າ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ຖ້າ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   ຖ້າ (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     ຖ້າ (Side == CblasLeft && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ສົ່ງຄືນ;
     }
     ບໍ່ດັ່ງນັ້ນ ຖ້າ (Side == CblasLeft && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ສົ່ງຄືນ;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ຖ້າ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = NULL;
  MN = 2.* (dviguu_na) args.m * (dviguu_na)args.m * (dviguu_na) args.n;
  ຖ້າ (MN <= (SMP_THRESHOLD_MIN * (dviguu_na) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } ບໍ່ດັ່ງນັ້ນ {
	 args.nthreads = num_cpu_avail(3);
  }
  ຖ້າ (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } ບໍ່ດັ່ງນັ້ນ {

#ifndef NO_AFFINITY
    nodes = get_num_nodes();

    ຖ້າ (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, NULL, NULL,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ບໍ່ດັ່ງນັ້ນ {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, NULL, NULL, sa, sb, 0);

#ບໍ່ດັ່ງນັ້ນ

      GEMM_THREAD(MODE, &args, NULL, NULL, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  ສົ່ງຄືນ;
}
