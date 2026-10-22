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

#شامل <stdio.h>
#شامل "common.h"
#شامل "l1param.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#وگرنه
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
ایستا صحیح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#وگرنهاگر defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#وگرنه
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 بیرونی نویسه* gotoblas_corename(تهی);
#endif

#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ایستا درونخطی صحیح get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, صحیح ncpu) {
  #ifdef DOUBLE
      بازگشت (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #وگرنه
      بازگشت
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ایستا درونخطی صحیح get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, صحیح ncpu) {
  بازگشت
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

ایستا درونخطی صحیح get_gemv_optimal_nthreads(BLASLONG MN) {
  صحیح ncpu = رقم_cpu_avail(3);
#اگر defined(_WIN64) && defined(_M_ARM64)
  اگر (MN > 100000000L)
    بازگشت رقم_cpu_avail(4);
  بازگشت 1;
#endif
#اگر defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  بازگشت get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#وگرنهاگر defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  بازگشت get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#وگرنهاگر defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  اگر (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    بازگشت get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  اگر (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    بازگشت get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  اگر ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    بازگشت 1;
  وگرنه
    بازگشت رقم_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  نویسه trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  صحیح buffer_size;
#ifdef SMP
  صحیح nthreads;
#endif

  صحیح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  اگر (trans == 'N') i = 0;
  اگر (trans == 'T') i = 1;
  اگر (trans == 'R') i = 0;
  اگر (trans == 'C') i = 1;

  اگر (incy == 0)	info = 11;
  اگر (incx == 0)	info = 8;
  اگر (lda < MAX(1, m))	info = 6;
  اگر (n < 0)		info = 3;
  اگر (m < 0)		info = 2;
  اگر (i < 0)          info = 1;

  trans = i;

  اگر (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order,
	   شمارش CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  صحیح trans, buffer_size;
  blasint info, t;
#ifdef SMP
  صحیح nthreads;
#endif

  صحیح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  اگر (order == CblasColMajor) {
    اگر (TransA == CblasNoTrans)     trans = 0;
    اگر (TransA == CblasTrans)       trans = 1;
    اگر (TransA == CblasConjNoTrans) trans = 0;
    اگر (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    اگر (incy == 0)	  info = 11;
    اگر (incx == 0)	  info = 8;
    اگر (lda < MAX(1, m))  info = 6;
    اگر (n < 0)		  info = 3;
    اگر (m < 0)		  info = 2;
    اگر (trans < 0)        info = 1;

  }

  اگر (order == CblasRowMajor) {
    اگر (TransA == CblasNoTrans)     trans = 1;
    اگر (TransA == CblasTrans)       trans = 0;
    اگر (TransA == CblasConjNoTrans) trans = 1;
    اگر (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    اگر (incy == 0)	  info = 11;
    اگر (incx == 0)	  info = 8;
    اگر (lda < MAX(1, m))  info = 6;
    اگر (n < 0)		  info = 3;
    اگر (m < 0)		  info = 2;
    اگر (trans < 0)        info = 1;

  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#endif
  اگر ((m==0) || (n==0)) بازگشت;

  lenx = n;
  leny = m;
  اگر (trans) lenx = m;
  اگر (trans) leny = n;

  اگر (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), پوچ, 0, پوچ, 0);

  اگر (alpha == ZERO) بازگشت;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  اگر (incx < 0) x -= (lenx - 1) * incx;
  اگر (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / اندازه(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / اندازه(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  پشته_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  اگر (nthreads == 1) {
#endif

    (gemv[(صحیح)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } وگرنه {

    (gemv_thread[(صحیح)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  پشته_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  بازگشت;

}
