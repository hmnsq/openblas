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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#அந்தர்பாவ "l1param.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#இல்லையெனில்
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
நிலையான முழு (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, முழு) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#அல்லதெனில் defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#இல்லையெனில்
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 வெளி எழுத்து* gotoblas_corename(வெறுமை);
#endif

#எனில் defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
நிலையான உள்வரி முழு get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, முழு ncpu) {
  #ifdef DOUBLE
      திரும்பு (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #இல்லையெனில்
      திரும்பு
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#எனில் defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
நிலையான உள்வரி முழு get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, முழு ncpu) {
  திரும்பு
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

நிலையான உள்வரி முழு get_gemv_optimal_nthreads(BLASLONG MN) {
  முழு ncpu = இலக்கம்_cpu_avail(3);
#எனில் defined(_WIN64) && defined(_M_ARM64)
  எனில் (MN > 100000000L)
    திரும்பு இலக்கம்_cpu_avail(4);
  திரும்பு 1;
#endif
#எனில் defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  திரும்பு get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#அல்லதெனில் defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  திரும்பு get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#அல்லதெனில் defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  எனில் (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    திரும்பு get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  எனில் (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    திரும்பு get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  எனில் ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    திரும்பு 1;
  இல்லையெனில்
    திரும்பு இலக்கம்_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  எழுத்து trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  முழு buffer_size;
#ifdef SMP
  முழு nthreads;
#endif

  முழு (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  எனில் (trans == 'N') i = 0;
  எனில் (trans == 'T') i = 1;
  எனில் (trans == 'R') i = 0;
  எனில் (trans == 'C') i = 1;

  எனில் (incy == 0)	info = 11;
  எனில் (incx == 0)	info = 8;
  எனில் (lda < MAX(1, m))	info = 6;
  எனில் (n < 0)		info = 3;
  எனில் (m < 0)		info = 2;
  எனில் (i < 0)          info = 1;

  trans = i;

  எனில் (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  முழு trans, buffer_size;
  blasint info, t;
#ifdef SMP
  முழு nthreads;
#endif

  முழு (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (TransA == CblasNoTrans)     trans = 0;
    எனில் (TransA == CblasTrans)       trans = 1;
    எனில் (TransA == CblasConjNoTrans) trans = 0;
    எனில் (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    எனில் (incy == 0)	  info = 11;
    எனில் (incx == 0)	  info = 8;
    எனில் (lda < MAX(1, m))  info = 6;
    எனில் (n < 0)		  info = 3;
    எனில் (m < 0)		  info = 2;
    எனில் (trans < 0)        info = 1;

  }

  எனில் (order == CblasRowMajor) {
    எனில் (TransA == CblasNoTrans)     trans = 1;
    எனில் (TransA == CblasTrans)       trans = 0;
    எனில் (TransA == CblasConjNoTrans) trans = 1;
    எனில் (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    எனில் (incy == 0)	  info = 11;
    எனில் (incx == 0)	  info = 8;
    எனில் (lda < MAX(1, m))  info = 6;
    எனில் (n < 0)		  info = 3;
    எனில் (m < 0)		  info = 2;
    எனில் (trans < 0)        info = 1;

  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif
  எனில் ((m==0) || (n==0)) திரும்பு;

  lenx = n;
  leny = m;
  எனில் (trans) lenx = m;
  எனில் (trans) leny = n;

  எனில் (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), வெற்று, 0, வெற்று, 0);

  எனில் (alpha == ZERO) திரும்பு;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  எனில் (incx < 0) x -= (lenx - 1) * incx;
  எனில் (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / அளவு(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / அளவு(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  அடுக்கு_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  எனில் (nthreads == 1) {
#endif

    (gemv[(முழு)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } இல்லையெனில் {

    (gemv_thread[(முழு)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  அடுக்கு_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  திரும்பு;

}
