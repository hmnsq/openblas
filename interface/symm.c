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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <ctype.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#નહીંતર
#define ERROR_NAME "SSYMM "
#endif
#નહીંતર
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#નહીંતર
#define ERROR_NAME "CSYMM "
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#નહીંતર
#define ERROR_NAME "CHEMM "
#endif
#endif
#નહીંતર
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#નહીંતર
#define ERROR_NAME "CSYMM3M "
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#નહીંતર
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
બાહ્ય અક્ષર* gotoblas_corename(શૂન્ય);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#અથવાજો defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#નહીંતર
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#નહીંતર
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#અથવાજો defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#નહીંતર
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

સ્થિર પૂર્ણાંક (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#નહીંતર
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#નહીંતર
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#નહીંતર
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *SIDE, અક્ષર *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  અક્ષર side_arg  = *SIDE;
  અક્ષર uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#જો defined(SMP) && !defined(NO_AFFINITY)
  પૂર્ણાંક nodes;
#endif
# જો defined(SMP)
  દ્વિગુણ MN;
#endif
  blasint info;
  પૂર્ણાંક side;
  પૂર્ણાંક uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (શૂન્ય *)alpha;
  args.beta  = (શૂન્ય *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  જો (side_arg  == 'L') side  = 0;
  જો (side_arg  == 'R') side  = 1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (શૂન્ય *)c;
  args.ldc = *ldC;

  info = 0;

  જો (args.ldc < MAX(1, args.m)) info = 12;

  જો (!side) {
    args.a = (શૂન્ય *)a;
    args.b = (શૂન્ય *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    જો (args.ldb < MAX(1, args.m)) info =  9;
    જો (args.lda < MAX(1, args.m)) info =  7;

  } નહીંતર {
    args.a = (શૂન્ય *)b;
    args.b = (શૂન્ય *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  જો (args.lda < MAX(1, args.m)) info =  9;
  જો (args.ldb < MAX(1, args.n)) info =  7;
  }

  જો (args.n   < 0)              info =  4;
  જો (args.m   < 0)              info =  3;
  જો (uplo     < 0)              info =  2;
  જો (side     < 0)              info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_SIDE Side, ગણના CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#નહીંતર
	   શૂન્ય *valpha,
	   શૂન્ય *va, blasint lda,
	   શૂન્ય *vb, blasint ldb,
	   શૂન્ય *vbeta,
	   શૂન્ય *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  પૂર્ણાંક side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#જો defined(SMP) && !defined(NO_AFFINITY)
  પૂર્ણાંક nodes;
#endif
#જો defined(SMP)
  દ્વિગુણ MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (શૂન્ય *)&alpha;
  args.beta  = (શૂન્ય *)&beta;
#નહીંતર
  args.alpha = (શૂન્ય *)alpha;
  args.beta  = (શૂન્ય *)beta;
#endif

  args.c = (શૂન્ય *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (Side == Cblas_ડાબું)  side = 0;
    જો (Side == Cblas_જમણું) side = 1;

    જો (Uplo == CblasUpper) uplo  = 0;
    જો (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    જો (args.ldc < MAX(1, args.m)) info = 12;

    જો (!side) {
      args.a = (શૂન્ય *)a;
      args.b = (શૂન્ય *)b;

      args.lda = lda;
      args.ldb = ldb;

      જો (args.ldb < MAX(1, args.m)) info =  9;
      જો (args.lda < MAX(1, args.m)) info =  7;

    } નહીંતર {
      args.a = (શૂન્ય *)b;
      args.b = (શૂન્ય *)a;

      args.lda = ldb;
      args.ldb = lda;

      જો (args.lda < MAX(1, args.m)) info =  9;
      જો (args.ldb < MAX(1, args.n)) info =  7;
    }

    જો (args.n   < 0)              info =  4;
    જો (args.m   < 0)              info =  3;
    જો (uplo     < 0)              info =  2;
    જો (side     < 0)              info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Side == Cblas_ડાબું)  side = 1;
    જો (Side == Cblas_જમણું) side = 0;

    જો (Uplo == CblasUpper) uplo  = 1;
    જો (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    જો (args.ldc < MAX(1, args.m)) info = 12;

    જો (!side) {
      args.a = (શૂન્ય *)a;
      args.b = (શૂન્ય *)b;

      args.lda = lda;
      args.ldb = ldb;

      જો (args.ldb < MAX(1, args.m)) info =  9;
      જો (args.lda < MAX(1, args.m)) info =  7;

    } નહીંતર {
      args.a = (શૂન્ય *)b;
      args.b = (શૂન્ય *)a;

      args.lda = ldb;
      args.ldb = lda;

      જો (args.lda < MAX(1, args.m)) info =  9;
      જો (args.ldb < MAX(1, args.n)) info =  7;
    }

    જો (args.n   < 0)              info =  4;
    જો (args.m   < 0)              info =  3;
    જો (uplo     < 0)              info =  2;
    જો (side     < 0)              info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

   જો (args.m == 0 || args.n == 0) પાછા;
#જો !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#જો defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#જો defined(DYNAMIC_ARCH)
જો (strcmp(gotoblas_corename(), "armv9sme") == 0
#જો defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   જો (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     જો (Side == Cblas_ડાબું && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); પાછા;
     }
     નહીંતર જો (Side == Cblas_ડાબું && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); પાછા;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = ખાલી;
  MN = 2.* (દ્વિગુણ) args.m * (દ્વિગુણ)args.m * (દ્વિગુણ) args.n;
  જો (MN <= (SMP_THRESHOLD_MIN * (દ્વિગુણ) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } નહીંતર {
	 args.nthreads = અંક_cpu_avail(3);
  }
  જો (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, ખાલી, ખાલી, sa, sb, 0);

#ifdef SMP
  } નહીંતર {

#ifndef NO_AFFINITY
    nodes = get_અંક_nodes();

    જો (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, ખાલી, ખાલી,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } નહીંતર {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, ખાલી, ખાલી, sa, sb, 0);

#નહીંતર

      GEMM_THREAD(MODE, &args, ખાલી, ખાલી, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  પાછા;
}
