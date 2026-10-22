/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project                               */
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
#include "common.h"
#include "l1param.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#иначе
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
static целое (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, целое) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#elif defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#иначе
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 extern char* gotoblas_corename(void);
#endif

#если defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
static inline целое get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, целое ncpu) {
  #ifdef DOUBLE
      вернуть (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #иначе
      вернуть
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#если defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
static inline целое get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, целое ncpu) {
  вернуть
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

static inline целое get_gemv_optimal_nthreads(BLASLONG MN) {
  целое ncpu = num_cpu_avail(3);
#если defined(_WIN64) && defined(_M_ARM64)
  если (MN > 100000000L)
    вернуть num_cpu_avail(4);
  вернуть 1;
#endif
#если defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  вернуть get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#elif defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  вернуть get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#elif defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  если (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    вернуть get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  если (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    вернуть get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  если ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    вернуть 1;
  иначе
    вернуть num_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  char trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  целое buffer_size;
#ifdef SMP
  целое nthreads;
#endif

  целое (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  если (trans == 'N') i = 0;
  если (trans == 'T') i = 1;
  если (trans == 'R') i = 0;
  если (trans == 'C') i = 1;

  если (incy == 0)	info = 11;
  если (incx == 0)	info = 8;
  если (lda < MAX(1, m))	info = 6;
  если (n < 0)		info = 3;
  если (m < 0)		info = 2;
  если (i < 0)          info = 1;

  trans = i;

  если (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#иначе

void CNAME(enum CBLAS_ORDER order,
	   enum CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  целое trans, buffer_size;
  blasint info, t;
#ifdef SMP
  целое nthreads;
#endif

  целое (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  если (order == CblasColMajor) {
    если (TransA == CblasNoTrans)     trans = 0;
    если (TransA == CblasTrans)       trans = 1;
    если (TransA == CblasConjNoTrans) trans = 0;
    если (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    если (incy == 0)	  info = 11;
    если (incx == 0)	  info = 8;
    если (lda < MAX(1, m))  info = 6;
    если (n < 0)		  info = 3;
    если (m < 0)		  info = 2;
    если (trans < 0)        info = 1;

  }

  если (order == CblasRowMajor) {
    если (TransA == CblasNoTrans)     trans = 1;
    если (TransA == CblasTrans)       trans = 0;
    если (TransA == CblasConjNoTrans) trans = 1;
    если (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    если (incy == 0)	  info = 11;
    если (incx == 0)	  info = 8;
    если (lda < MAX(1, m))  info = 6;
    если (n < 0)		  info = 3;
    если (m < 0)		  info = 2;
    если (trans < 0)        info = 1;

  }

  если (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    вернуть;
  }

#endif
  если ((m==0) || (n==0)) вернуть;

  lenx = n;
  leny = m;
  если (trans) lenx = m;
  если (trans) leny = n;

  если (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  если (alpha == ZERO) вернуть;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  если (incx < 0) x -= (lenx - 1) * incx;
  если (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / sizeof(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / sizeof(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  если (nthreads == 1) {
#endif

    (gemv[(целое)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } иначе {

    (gemv_thread[(целое)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  вернуть;

}
