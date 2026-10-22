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
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DSYMM "
#अञ्ञथा
#define ERROR_NAME "SSYMM "
#endif
#अञ्ञथा
#define SMP_THRESHOLD_MIN 8192.
#ifndef GEMM3M
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYMM "
#अञ्ञथा
#define ERROR_NAME "CSYMM "
#endif
#अञ्ञथा
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHEMM "
#अञ्ञथा
#define ERROR_NAME "CHEMM "
#endif
#endif
#अञ्ञथा
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYMM3M "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYMM3M "
#अञ्ञथा
#define ERROR_NAME "CSYMM3M "
#endif
#अञ्ञथा
#ifdef XDOUBLE
#define ERROR_NAME "XHEMM3M "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHEMM3M "
#अञ्ञथा
#define ERROR_NAME "CHEMM3M "
#endif
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

#ifdef DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);                                
#endif                                                               

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_REAL)
#नोचेद्यदि defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_REAL)
#अञ्ञथा
#define MODE	(BLAS_SINGLE  | BLAS_REAL)
#endif
#अञ्ञथा
#ifdef XDOUBLE
#define MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#नोचेद्यदि defined(DOUBLE)
#define MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#अञ्ञथा
#define MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#endif
#endif
#endif

स्थिर पूर्णाङ्क (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef GEMM3M
#ifndef HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#सचे defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#endif
#अञ्ञथा
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#सचे defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#endif
#endif
#अञ्ञथा
#ifndef HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#सचे defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#endif
#अञ्ञथा
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#सचे defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#endif
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *SIDE, वर्ण *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण side_arg  = *SIDE;
  वर्ण uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#सचे defined(SMP) && !defined(NO_AFFINITY)
  पूर्णाङ्क nodes;
#endif
# सचे defined(SMP)
  द्विगुण MN;
#endif
  blasint info;
  पूर्णाङ्क side;
  पूर्णाङ्क uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  सचे (side_arg  == 'L') side  = 0;
  सचे (side_arg  == 'R') side  = 1;

  सचे (uplo_arg  == 'U') uplo  = 0;
  सचे (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (शून्य *)c;
  args.ldc = *ldC;

  info = 0;

  सचे (args.ldc < MAX(1, args.m)) info = 12;

  सचे (!side) {
    args.a = (शून्य *)a;
    args.b = (शून्य *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    सचे (args.ldb < MAX(1, args.m)) info =  9;
    सचे (args.lda < MAX(1, args.m)) info =  7;

  } अञ्ञथा {
    args.a = (शून्य *)b;
    args.b = (शून्य *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  सचे (args.lda < MAX(1, args.m)) info =  9;
  सचे (args.ldb < MAX(1, args.n)) info =  7;
  }

  सचे (args.n   < 0)              info =  4;
  सचे (args.m   < 0)              info =  3;
  सचे (uplo     < 0)              info =  2;
  सचे (side     < 0)              info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_SIDE Side, गणना CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अञ्ञथा
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
  पूर्णाङ्क side, uplo;
  blasint info;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#सचे defined(SMP) && !defined(NO_AFFINITY)
  पूर्णाङ्क nodes;
#endif
#सचे defined(SMP)
  द्विगुण MN;
#endif

  PRINT_DEBUG_CNAME;

#ifndef COMPLEX
  args.alpha = (शून्य *)&alpha;
  args.beta  = (शून्य *)&beta;
#अञ्ञथा
  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;
#endif

  args.c = (शून्य *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  सचे (order == CblasColMajor) {
    सचे (Side == Cblas_वामः)  side = 0;
    सचे (Side == Cblas_दक्षिणः) side = 1;

    सचे (Uplo == CblasUpper) uplo  = 0;
    सचे (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    सचे (args.ldc < MAX(1, args.m)) info = 12;

    सचे (!side) {
      args.a = (शून्य *)a;
      args.b = (शून्य *)b;

      args.lda = lda;
      args.ldb = ldb;

      सचे (args.ldb < MAX(1, args.m)) info =  9;
      सचे (args.lda < MAX(1, args.m)) info =  7;

    } अञ्ञथा {
      args.a = (शून्य *)b;
      args.b = (शून्य *)a;

      args.lda = ldb;
      args.ldb = lda;

      सचे (args.lda < MAX(1, args.m)) info =  9;
      सचे (args.ldb < MAX(1, args.n)) info =  7;
    }

    सचे (args.n   < 0)              info =  4;
    सचे (args.m   < 0)              info =  3;
    सचे (uplo     < 0)              info =  2;
    सचे (side     < 0)              info =  1;
  }

  सचे (order == CblasRowMajor) {
    सचे (Side == Cblas_वामः)  side = 1;
    सचे (Side == Cblas_दक्षिणः) side = 0;

    सचे (Uplo == CblasUpper) uplo  = 1;
    सचे (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    सचे (args.ldc < MAX(1, args.m)) info = 12;

    सचे (!side) {
      args.a = (शून्य *)a;
      args.b = (शून्य *)b;

      args.lda = lda;
      args.ldb = ldb;

      सचे (args.ldb < MAX(1, args.m)) info =  9;
      सचे (args.lda < MAX(1, args.m)) info =  7;

    } अञ्ञथा {
      args.a = (शून्य *)b;
      args.b = (शून्य *)a;

      args.lda = ldb;
      args.ldb = lda;

      सचे (args.lda < MAX(1, args.m)) info =  9;
      सचे (args.ldb < MAX(1, args.n)) info =  7;
    }

    सचे (args.n   < 0)              info =  4;
    सचे (args.m   < 0)              info =  3;
    सचे (uplo     < 0)              info =  2;
    सचे (side     < 0)              info =  1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

   सचे (args.m == 0 || args.n == 0) पटिनिवत्त;
#सचे !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#सचे defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#सचे defined(DYNAMIC_ARCH)
सचे (strcmp(gotoblas_corename(), "armv9sme") == 0
#सचे defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   सचे (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     सचे (Side == Cblas_वामः && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); पटिनिवत्त;
     }
     अञ्ञथा सचे (Side == Cblas_वामः && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); पटिनिवत्त;
     }
   }
#endif
#endif

#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  सचे (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  args.common = रिक्त;
  MN = 2.* (द्विगुण) args.m * (द्विगुण)args.m * (द्विगुण) args.n;
  सचे (MN <= (SMP_THRESHOLD_MIN * (द्विगुण) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } अञ्ञथा {
	 args.nthreads = अङ्कः_cpu_avail(3);
  }
  सचे (args.nthreads == 1) {
#endif

    (symm[(side << 1) | uplo ])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP
  } अञ्ञथा {

#ifndef NO_AFFINITY
    nodes = get_अङ्कः_nodes();

    सचे (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, रिक्त, रिक्त,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } अञ्ञथा {
#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, रिक्त, रिक्त, sa, sb, 0);

#अञ्ञथा

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

  पटिनिवत्त;
}
