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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ଅନ୍ତର୍ଭାଵ "l1param.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#ନଚେତଯଦି defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ଅନ୍ୟଥା
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 ବାହ୍ୟ ଅକ୍ଷର* gotoblas_corename(ଶୂନ୍ୟ);
#endif

#ଯଦି defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu) {
  #ifdef DOUBLE
      ଫେରନ୍ତୁ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ଅନ୍ୟଥା
      ଫେରନ୍ତୁ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#ଯଦି defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu) {
  ଫେରନ୍ତୁ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

ସ୍ଥିର ଧାଡ଼ିଭିତରେ ପୂର୍ଣ୍ଣସଂଖ୍ୟା get_gemv_optimal_nthreads(BLASLONG MN) {
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା ncpu = ଅଙ୍କ_cpu_avail(3);
#ଯଦି defined(_WIN64) && defined(_M_ARM64)
  ଯଦି (MN > 100000000L)
    ଫେରନ୍ତୁ ଅଙ୍କ_cpu_avail(4);
  ଫେରନ୍ତୁ 1;
#endif
#ଯଦି defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ଫେରନ୍ତୁ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#ନଚେତଯଦି defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ଫେରନ୍ତୁ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#ନଚେତଯଦି defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ଯଦି (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ଫେରନ୍ତୁ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ଯଦି (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ଫେରନ୍ତୁ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  ଯଦି ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    ଫେରନ୍ତୁ 1;
  ଅନ୍ୟଥା
    ଫେରନ୍ତୁ ଅଙ୍କ_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  ଅକ୍ଷର trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା buffer_size;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ଯଦି (trans == 'N') i = 0;
  ଯଦି (trans == 'T') i = 1;
  ଯଦି (trans == 'R') i = 0;
  ଯଦି (trans == 'C') i = 1;

  ଯଦି (incy == 0)	info = 11;
  ଯଦି (incx == 0)	info = 8;
  ଯଦି (lda < MAX(1, m))	info = 6;
  ଯଦି (n < 0)		info = 3;
  ଯଦି (m < 0)		info = 2;
  ଯଦି (i < 0)          info = 1;

  trans = i;

  ଯଦି (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order,
	   ଗଣନା CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans, buffer_size;
  blasint info, t;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 0;
    ଯଦି (TransA == CblasTrans)       trans = 1;
    ଯଦି (TransA == CblasConjNoTrans) trans = 0;
    ଯଦି (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ଯଦି (incy == 0)	  info = 11;
    ଯଦି (incx == 0)	  info = 8;
    ଯଦି (lda < MAX(1, m))  info = 6;
    ଯଦି (n < 0)		  info = 3;
    ଯଦି (m < 0)		  info = 2;
    ଯଦି (trans < 0)        info = 1;

  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 1;
    ଯଦି (TransA == CblasTrans)       trans = 0;
    ଯଦି (TransA == CblasConjNoTrans) trans = 1;
    ଯଦି (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    ଯଦି (incy == 0)	  info = 11;
    ଯଦି (incx == 0)	  info = 8;
    ଯଦି (lda < MAX(1, m))  info = 6;
    ଯଦି (n < 0)		  info = 3;
    ଯଦି (m < 0)		  info = 2;
    ଯଦି (trans < 0)        info = 1;

  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif
  ଯଦି ((m==0) || (n==0)) ଫେରନ୍ତୁ;

  lenx = n;
  leny = m;
  ଯଦି (trans) lenx = m;
  ଯଦି (trans) leny = n;

  ଯଦି (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ଖାଲି, 0, ଖାଲି, 0);

  ଯଦି (alpha == ZERO) ଫେରନ୍ତୁ;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0) x -= (lenx - 1) * incx;
  ଯଦି (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / ଆକାର(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / ଆକାର(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ଗଦା_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ଯଦି (nthreads == 1) {
#endif

    (gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ଅନ୍ୟଥା {

    (gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  ଗଦା_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;

}
