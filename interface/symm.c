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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#অন্যথা
#define ERROR_NAME "SSYMM "
#endif
#অন্যথা
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#অন্যথা
#define ERROR_NAME "CSYMM "
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#অন্যথা
#define ERROR_NAME "CHEMM "
#endif
#endif
#অন্যথা
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#অন্যথা
#define ERROR_NAME "CSYMM3M "
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#অন্যথা
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
বাহ্যিক অক্ষর* gotoblas_corename(শূন্য);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#নতুবা_যদি defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#অন্যথা
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#অন্যথা
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#নতুবা_যদি defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#অন্যথা
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

স্থির পূর্ণসংখ্যা (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#অন্যথা
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#অন্যথা
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#অন্যথা
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(অক্ষর *SIDE, অক্ষর *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  অক্ষর side_arg  = *SIDE;
  অক্ষর uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *বাফার;
  FLOAT *sa, *sb;

#যদি defined(SMP) && !defined(NO_AFFINITY)
  পূর্ণসংখ্যা nodes;
#endif
# যদি defined(SMP)
  দ্বিগুণ MN;
#endif
  blasint info;
  পূর্ণসংখ্যা side;
  পূর্ণসংখ্যা uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  যদি (side_arg  == 'L') side  = 0;
  যদি (side_arg  == 'R') side  = 1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (শূন্য *)c;
  args.ldc = *ldC;

  info = 0;

  যদি (args.ldc < MAX(1, args.m)) info = 12;

  যদি (!side) {
    args.a = (শূন্য *)a;
    args.b = (শূন্য *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    যদি (args.ldb < MAX(1, args.m)) info =  9;
    যদি (args.lda < MAX(1, args.m)) info =  7;

  } অন্যথা {
    args.a = (শূন্য *)b;
    args.b = (শূন্য *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  যদি (args.lda < MAX(1, args.m)) info =  9;
  যদি (args.ldb < MAX(1, args.n)) info =  7;
  }

  যদি (args.n   < 0)              info =  4;
  যদি (args.m   < 0)              info =  3;
  যদি (uplo     < 0)              info =  2;
  যদি (side     < 0)              info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_SIDE Side, গণনা CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#অন্যথা
	   শূন্য *valpha,
	   শূন্য *va, blasint lda,
	   শূন্য *vb, blasint ldb,
	   শূন্য *vbeta,
	   শূন্য *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  পূর্ণসংখ্যা side, uplo;
  blasint info;

  FLOAT *বাফার;
  FLOAT *sa, *sb;

#যদি defined(SMP) && !defined(NO_AFFINITY)
  পূর্ণসংখ্যা nodes;
#endif
#যদি defined(SMP)
  দ্বিগুণ MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (শূন্য *)&alpha;
  args.beta  = (শূন্য *)&beta;
#অন্যথা
  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;
#endif

  args.c = (শূন্য *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (Side == Cblas_বাম)  side = 0;
    যদি (Side == Cblas_ডান) side = 1;

    যদি (Uplo == CblasUpper) uplo  = 0;
    যদি (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    যদি (args.ldc < MAX(1, args.m)) info = 12;

    যদি (!side) {
      args.a = (শূন্য *)a;
      args.b = (শূন্য *)b;

      args.lda = lda;
      args.ldb = ldb;

      যদি (args.ldb < MAX(1, args.m)) info =  9;
      যদি (args.lda < MAX(1, args.m)) info =  7;

    } অন্যথা {
      args.a = (শূন্য *)b;
      args.b = (শূন্য *)a;

      args.lda = ldb;
      args.ldb = lda;

      যদি (args.lda < MAX(1, args.m)) info =  9;
      যদি (args.ldb < MAX(1, args.n)) info =  7;
    }

    যদি (args.n   < 0)              info =  4;
    যদি (args.m   < 0)              info =  3;
    যদি (uplo     < 0)              info =  2;
    যদি (side     < 0)              info =  1;
  }

  যদি (order == CblasRowMajor) {
    যদি (Side == Cblas_বাম)  side = 1;
    যদি (Side == Cblas_ডান) side = 0;

    যদি (Uplo == CblasUpper) uplo  = 1;
    যদি (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    যদি (args.ldc < MAX(1, args.m)) info = 12;

    যদি (!side) {
      args.a = (শূন্য *)a;
      args.b = (শূন্য *)b;

      args.lda = lda;
      args.ldb = ldb;

      যদি (args.ldb < MAX(1, args.m)) info =  9;
      যদি (args.lda < MAX(1, args.m)) info =  7;

    } অন্যথা {
      args.a = (শূন্য *)b;
      args.b = (শূন্য *)a;

      args.lda = ldb;
      args.ldb = lda;

      যদি (args.lda < MAX(1, args.m)) info =  9;
      যদি (args.ldb < MAX(1, args.n)) info =  7;
    }

    যদি (args.n   < 0)              info =  4;
    যদি (args.m   < 0)              info =  3;
    যদি (uplo     < 0)              info =  2;
    যদি (side     < 0)              info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

   যদি (args.m == 0 || args.n == 0) ফেরত;
#যদি !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#যদি defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
যদি (strcmp(gotoblas_corename(), "armv9sme") == 0
#যদি defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   যদি (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     যদি (Side == Cblas_বাম && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ফেরত;
     }
     অন্যথা যদি (Side == Cblas_বাম && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); ফেরত;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  বাফার = (FLOAT *)blas_memory_alloc(0);
  যদি (!বাফার) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

  sa = (FLOAT *)((BLASLONG)বাফার + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * আকার_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = ফাঁকা;
  MN = 2.* (দ্বিগুণ) args.m * (দ্বিগুণ)args.m * (দ্বিগুণ) args.n;
  যদি (MN <= (SMP_THRESHOLD_MIN * (দ্বিগুণ) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } অন্যথা {
	 args.nthreads = অঙ্ক_cpu_avail(3);
  }
  যদি (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, ফাঁকা, ফাঁকা, sa, sb, 0);

#ifdef SMP
  } অন্যথা {

#ifndef NO_AFFINITY
    nodes = get_অঙ্ক_nodes();

    যদি (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, ফাঁকা, ফাঁকা,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } অন্যথা {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, ফাঁকা, ফাঁকা, sa, sb, 0);

#অন্যথা

      GEMM_THREAD(MODE, &args, ফাঁকা, ফাঁকা, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

#endif

#ifndef NO_AFFINITY
    }
#endif

  }
#endif

 blas_memory_free(বাফার);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side)? args.m * (args.m / 2 + args.n) : args.n * (args.m + args.n / 2),
		       (!side)? 2 * args.m * args.m * args.n : 2 * args.m * args.n * args.n);

  IDEBUG_END;

  ফেরত;
}
