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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#অন্তর্ভাৱ "l1param.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#অন্যথা
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
স্থির পূর্ণাঙ্ক (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#নোচেদ্যদি defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#অন্যথা
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 বাহ্য ৱর্ণ* gotoblas_corename(শূন্য);
#endif

#যদি defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, পূর্ণাঙ্ক ncpu) {
  #ifdef DOUBLE
      ঘূৰাই (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #অন্যথা
      ঘূৰাই
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#যদি defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, পূর্ণাঙ্ক ncpu) {
  ঘূৰাই
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

স্থির অন্তঃস্থ পূর্ণাঙ্ক get_gemv_optimal_nthreads(BLASLONG MN) {
  পূর্ণাঙ্ক ncpu = অঙ্ক_cpu_avail(3);
#যদি defined(_WIN64) && defined(_M_ARM64)
  যদি (MN > 100000000L)
    ঘূৰাই অঙ্ক_cpu_avail(4);
  ঘূৰাই 1;
#endif
#যদি defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ঘূৰাই get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#নোচেদ্যদি defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ঘূৰাই get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#নোচেদ্যদি defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  যদি (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ঘূৰাই get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  যদি (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ঘূৰাই get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  যদি ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    ঘূৰাই 1;
  অন্যথা
    ঘূৰাই অঙ্ক_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(ৱর্ণ *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  ৱর্ণ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  পূর্ণাঙ্ক buffer_size;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  পূর্ণাঙ্ক (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  যদি (trans == 'N') i = 0;
  যদি (trans == 'T') i = 1;
  যদি (trans == 'R') i = 0;
  যদি (trans == 'C') i = 1;

  যদি (incy == 0)	info = 11;
  যদি (incx == 0)	info = 8;
  যদি (lda < MAX(1, m))	info = 6;
  যদি (n < 0)		info = 3;
  যদি (m < 0)		info = 2;
  যদি (i < 0)          info = 1;

  trans = i;

  যদি (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  পূর্ণাঙ্ক trans, buffer_size;
  blasint info, t;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  পূর্ণাঙ্ক (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 0;
    যদি (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    যদি (incy == 0)	  info = 11;
    যদি (incx == 0)	  info = 8;
    যদি (lda < MAX(1, m))  info = 6;
    যদি (n < 0)		  info = 3;
    যদি (m < 0)		  info = 2;
    যদি (trans < 0)        info = 1;

  }

  যদি (order == CblasRowMajor) {
    যদি (TransA == CblasNoTrans)     trans = 1;
    যদি (TransA == CblasTrans)       trans = 0;
    যদি (TransA == CblasConjNoTrans) trans = 1;
    যদি (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    যদি (incy == 0)	  info = 11;
    যদি (incx == 0)	  info = 8;
    যদি (lda < MAX(1, m))  info = 6;
    যদি (n < 0)		  info = 3;
    যদি (m < 0)		  info = 2;
    যদি (trans < 0)        info = 1;

  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    ঘূৰাই;
  }

#endif
  যদি ((m==0) || (n==0)) ঘূৰাই;

  lenx = n;
  leny = m;
  যদি (trans) lenx = m;
  যদি (trans) leny = n;

  যদি (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);

  যদি (alpha == ZERO) ঘূৰাই;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (lenx - 1) * incx;
  যদি (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / পরিমাণ(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / পরিমাণ(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  স্তূপ_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  যদি (nthreads == 1) {
#endif

    (gemv[(পূর্ণাঙ্ক)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } অন্যথা {

    (gemv_thread[(পূর্ণাঙ্ক)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  স্তূপ_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ঘূৰাই;

}
