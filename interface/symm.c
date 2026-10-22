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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.
#ifdef XDOUBLE
#define ERROR_NAME "QSYMM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#अन्यथा
#define ERROR_NAME "SSYMM "
#endif
#अन्यथा
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#अन्यथा
#define ERROR_NAME "CSYMM "
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#अन्यथा
#define ERROR_NAME "CHEMM "
#endif
#endif
#अन्यथा
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#अन्यथा
#define ERROR_NAME "CSYMM3M "
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#अन्यथा
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(शून्य);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#किंवाजर defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#अन्यथा
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#अन्यथा
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#किंवाजर defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#अन्यथा
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

स्थिर पूर्णांक (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#अन्यथा
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#अन्यथा
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#अन्यथा
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#जर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *SIDE, अक्षर *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  अक्षर side_arg  = *SIDE;
  अक्षर uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#जर defined(SMP) && !defined(NO_AFFINITY)
  पूर्णांक nodes;
#endif
# जर defined(SMP)
  द्विगुण MN;
#endif
  blasint info;
  पूर्णांक side;
  पूर्णांक uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  जर (side_arg  == 'L') side  = 0;
  जर (side_arg  == 'R') side  = 1;

  जर (uplo_arg  == 'U') uplo  = 0;
  जर (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (शून्य *)c;
  args.ldc = *ldC;

  info = 0;

  जर (args.ldc < MAX(1, args.m)) info = 12;

  जर (!side) {
    args.a = (शून्य *)a;
    args.b = (शून्य *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    जर (args.ldb < MAX(1, args.m)) info =  9;
    जर (args.lda < MAX(1, args.m)) info =  7;

  } अन्यथा {
    args.a = (शून्य *)b;
    args.b = (शून्य *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  जर (args.lda < MAX(1, args.m)) info =  9;
  जर (args.ldb < MAX(1, args.n)) info =  7;
  }

  जर (args.n   < 0)              info =  4;
  जर (args.m   < 0)              info =  3;
  जर (uplo     < 0)              info =  2;
  जर (side     < 0)              info =  1;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_SIDE Side, गणना CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
	   शून्य *vbeta,
	   शून्य *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#endif

  blas_arg_t args;
  पूर्णांक side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#जर defined(SMP) && !defined(NO_AFFINITY)
  पूर्णांक nodes;
#endif
#जर defined(SMP)
  द्विगुण MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (शून्य *)&alpha;
  args.beta  = (शून्य *)&beta;
#अन्यथा
  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;
#endif

  args.c = (शून्य *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  जर (order == CblasColMajor) {
    जर (Side == Cblas_डावा)  side = 0;
    जर (Side == Cblas_उजवा) side = 1;

    जर (Uplo == CblasUpper) uplo  = 0;
    जर (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    जर (args.ldc < MAX(1, args.m)) info = 12;

    जर (!side) {
      args.a = (शून्य *)a;
      args.b = (शून्य *)b;

      args.lda = lda;
      args.ldb = ldb;

      जर (args.ldb < MAX(1, args.m)) info =  9;
      जर (args.lda < MAX(1, args.m)) info =  7;

    } अन्यथा {
      args.a = (शून्य *)b;
      args.b = (शून्य *)a;

      args.lda = ldb;
      args.ldb = lda;

      जर (args.lda < MAX(1, args.m)) info =  9;
      जर (args.ldb < MAX(1, args.n)) info =  7;
    }

    जर (args.n   < 0)              info =  4;
    जर (args.m   < 0)              info =  3;
    जर (uplo     < 0)              info =  2;
    जर (side     < 0)              info =  1;
  }

  जर (order == CblasRowMajor) {
    जर (Side == Cblas_डावा)  side = 1;
    जर (Side == Cblas_उजवा) side = 0;

    जर (Uplo == CblasUpper) uplo  = 1;
    जर (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    जर (args.ldc < MAX(1, args.m)) info = 12;

    जर (!side) {
      args.a = (शून्य *)a;
      args.b = (शून्य *)b;

      args.lda = lda;
      args.ldb = ldb;

      जर (args.ldb < MAX(1, args.m)) info =  9;
      जर (args.lda < MAX(1, args.m)) info =  7;

    } अन्यथा {
      args.a = (शून्य *)b;
      args.b = (शून्य *)a;

      args.lda = ldb;
      args.ldb = lda;

      जर (args.lda < MAX(1, args.m)) info =  9;
      जर (args.ldb < MAX(1, args.n)) info =  7;
    }

    जर (args.n   < 0)              info =  4;
    जर (args.m   < 0)              info =  3;
    जर (uplo     < 0)              info =  2;
    जर (side     < 0)              info =  1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

   जर (args.m == 0 || args.n == 0) परत;
#जर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#जर defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
जर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#जर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   जर (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     जर (Side == Cblas_डावा && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); परत;
     }
     अन्यथा जर (Side == Cblas_डावा && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); परत;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = रिक्त;
  MN = 2.* (द्विगुण) args.m * (द्विगुण)args.m * (द्विगुण) args.n;
  जर (MN <= (SMP_THRESHOLD_MIN * (द्विगुण) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } अन्यथा {
	 args.nthreads = अंक_cpu_avail(3);
  }
  जर (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP
  } अन्यथा {

#ifndef NO_AFFINITY
    nodes = get_अंक_nodes();

    जर (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, रिक्त, रिक्त,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } अन्यथा {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

      GEMM_THREAD(MODE, &args, रिक्त, रिक्त, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

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

  परत;
}
