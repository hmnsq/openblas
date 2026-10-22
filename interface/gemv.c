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

#ضمن <stdio.h>
#ضمن "common.h"
#ضمن "l1param.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#ܐܚܪܢܐ
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
ساكن صحيح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#وإلاإذا defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ܐܚܪܢܐ
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 خارجي حرف* gotoblas_corename(فراغ);
#endif

#ܐܢ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ساكن مضمن صحيح get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, صحيح ncpu) {
  #ifdef DOUBLE
      ܦܢܐ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ܐܚܪܢܐ
      ܦܢܐ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#ܐܢ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ساكن مضمن صحيح get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, صحيح ncpu) {
  ܦܢܐ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

ساكن مضمن صحيح get_gemv_optimal_nthreads(BLASLONG MN) {
  صحيح ncpu = num_cpu_avail(3);
#ܐܢ defined(_WIN64) && defined(_M_ARM64)
  ܐܢ (MN > 100000000L)
    ܦܢܐ num_cpu_avail(4);
  ܦܢܐ 1;
#endif
#ܐܢ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ܦܢܐ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#وإلاإذا defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ܦܢܐ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#وإلاإذا defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ܐܢ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ܦܢܐ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ܐܢ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ܦܢܐ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  ܐܢ ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    ܦܢܐ 1;
  ܐܚܪܢܐ
    ܦܢܐ num_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  حرف trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  صحيح buffer_size;
#ifdef SMP
  صحيح nthreads;
#endif

  صحيح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ܐܢ (trans == 'N') i = 0;
  ܐܢ (trans == 'T') i = 1;
  ܐܢ (trans == 'R') i = 0;
  ܐܢ (trans == 'C') i = 1;

  ܐܢ (incy == 0)	info = 11;
  ܐܢ (incx == 0)	info = 8;
  ܐܢ (lda < MAX(1, m))	info = 6;
  ܐܢ (n < 0)		info = 3;
  ܐܢ (m < 0)		info = 2;
  ܐܢ (i < 0)          info = 1;

  trans = i;

  ܐܢ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  صحيح trans, buffer_size;
  blasint info, t;
#ifdef SMP
  صحيح nthreads;
#endif

  صحيح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ܐܢ (order == CblasColMajor) {
    ܐܢ (TransA == CblasNoTrans)     trans = 0;
    ܐܢ (TransA == CblasTrans)       trans = 1;
    ܐܢ (TransA == CblasConjNoTrans) trans = 0;
    ܐܢ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ܐܢ (incy == 0)	  info = 11;
    ܐܢ (incx == 0)	  info = 8;
    ܐܢ (lda < MAX(1, m))  info = 6;
    ܐܢ (n < 0)		  info = 3;
    ܐܢ (m < 0)		  info = 2;
    ܐܢ (trans < 0)        info = 1;

  }

  ܐܢ (order == CblasRowMajor) {
    ܐܢ (TransA == CblasNoTrans)     trans = 1;
    ܐܢ (TransA == CblasTrans)       trans = 0;
    ܐܢ (TransA == CblasConjNoTrans) trans = 1;
    ܐܢ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    ܐܢ (incy == 0)	  info = 11;
    ܐܢ (incx == 0)	  info = 8;
    ܐܢ (lda < MAX(1, m))  info = 6;
    ܐܢ (n < 0)		  info = 3;
    ܐܢ (m < 0)		  info = 2;
    ܐܢ (trans < 0)        info = 1;

  }

  ܐܢ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#endif
  ܐܢ ((m==0) || (n==0)) ܦܢܐ;

  lenx = n;
  leny = m;
  ܐܢ (trans) lenx = m;
  ܐܢ (trans) leny = n;

  ܐܢ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  ܐܢ (alpha == ZERO) ܦܢܐ;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ܐܢ (incx < 0) x -= (lenx - 1) * incx;
  ܐܢ (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / حجم(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / حجم(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ܐܢ (nthreads == 1) {
#endif

    (gemv[(صحيح)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ܐܚܪܢܐ {

    (gemv_thread[(صحيح)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ܦܢܐ;

}
