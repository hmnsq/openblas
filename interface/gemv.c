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
#không_thì
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
static số_nguyên (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, số_nguyên) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#elif defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#không_thì
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 extern char* gotoblas_corename(void);
#endif

#nếu defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
static inline số_nguyên get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, số_nguyên ncpu) {
  #ifdef DOUBLE
      trả_về (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #không_thì
      trả_về
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#nếu defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
static inline số_nguyên get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, số_nguyên ncpu) {
  trả_về
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

static inline số_nguyên get_gemv_optimal_nthreads(BLASLONG MN) {
  số_nguyên ncpu = num_cpu_avail(3);
#nếu defined(_WIN64) && defined(_M_ARM64)
  nếu (MN > 100000000L)
    trả_về num_cpu_avail(4);
  trả_về 1;
#endif
#nếu defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  trả_về get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#elif defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  trả_về get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#elif defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  nếu (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    trả_về get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  nếu (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    trả_về get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  nếu ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    trả_về 1;
  không_thì
    trả_về num_cpu_avail(2);
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
  số_nguyên buffer_size;
#ifdef SMP
  số_nguyên nthreads;
#endif

  số_nguyên (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  nếu (trans == 'N') i = 0;
  nếu (trans == 'T') i = 1;
  nếu (trans == 'R') i = 0;
  nếu (trans == 'C') i = 1;

  nếu (incy == 0)	info = 11;
  nếu (incx == 0)	info = 8;
  nếu (lda < MAX(1, m))	info = 6;
  nếu (n < 0)		info = 3;
  nếu (m < 0)		info = 2;
  nếu (i < 0)          info = 1;

  trans = i;

  nếu (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

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
  số_nguyên trans, buffer_size;
  blasint info, t;
#ifdef SMP
  số_nguyên nthreads;
#endif

  số_nguyên (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  nếu (order == CblasColMajor) {
    nếu (TransA == CblasNoTrans)     trans = 0;
    nếu (TransA == CblasTrans)       trans = 1;
    nếu (TransA == CblasConjNoTrans) trans = 0;
    nếu (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    nếu (incy == 0)	  info = 11;
    nếu (incx == 0)	  info = 8;
    nếu (lda < MAX(1, m))  info = 6;
    nếu (n < 0)		  info = 3;
    nếu (m < 0)		  info = 2;
    nếu (trans < 0)        info = 1;

  }

  nếu (order == CblasRowMajor) {
    nếu (TransA == CblasNoTrans)     trans = 1;
    nếu (TransA == CblasTrans)       trans = 0;
    nếu (TransA == CblasConjNoTrans) trans = 1;
    nếu (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    nếu (incy == 0)	  info = 11;
    nếu (incx == 0)	  info = 8;
    nếu (lda < MAX(1, m))  info = 6;
    nếu (n < 0)		  info = 3;
    nếu (m < 0)		  info = 2;
    nếu (trans < 0)        info = 1;

  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif
  nếu ((m==0) || (n==0)) trả_về;

  lenx = n;
  leny = m;
  nếu (trans) lenx = m;
  nếu (trans) leny = n;

  nếu (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  nếu (alpha == ZERO) trả_về;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  nếu (incx < 0) x -= (lenx - 1) * incx;
  nếu (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / sizeof(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / sizeof(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  nếu (nthreads == 1) {
#endif

    (gemv[(số_nguyên)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } không_thì {

    (gemv_thread[(số_nguyên)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  trả_về;

}
