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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ಅನ್ತರ್ಭಾವ "l1param.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#ನೋಚೇದ್ಯದಿ defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ಅನ್ಯಥಾ
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 ಬಾಹ್ಯ ವರ್ಣ* gotoblas_corename(ಶೂನ್ಯ);
#endif

#ಯದಿ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, ಪೂರ್ಣಾಙ್ಕ ncpu) {
  #ifdef DOUBLE
      ಪ್ರತ್ಯಾಗಮ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ಅನ್ಯಥಾ
      ಪ್ರತ್ಯಾಗಮ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#ಯದಿ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, ಪೂರ್ಣಾಙ್ಕ ncpu) {
  ಪ್ರತ್ಯಾಗಮ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

ಸ್ಥಿರ ಅನ್ತಃಸ್ಥ ಪೂರ್ಣಾಙ್ಕ get_gemv_optimal_nthreads(BLASLONG MN) {
  ಪೂರ್ಣಾಙ್ಕ ncpu = ಅಂಕೆ_cpu_avail(3);
#ಯದಿ defined(_WIN64) && defined(_M_ARM64)
  ಯದಿ (MN > 100000000L)
    ಪ್ರತ್ಯಾಗಮ ಅಂಕೆ_cpu_avail(4);
  ಪ್ರತ್ಯಾಗಮ 1;
#endif
#ಯದಿ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ಪ್ರತ್ಯಾಗಮ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#ನೋಚೇದ್ಯದಿ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ಪ್ರತ್ಯಾಗಮ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#ನೋಚೇದ್ಯದಿ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ಯದಿ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ಪ್ರತ್ಯಾಗಮ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ಯದಿ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ಪ್ರತ್ಯಾಗಮ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  ಯದಿ ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    ಪ್ರತ್ಯಾಗಮ 1;
  ಅನ್ಯಥಾ
    ಪ್ರತ್ಯಾಗಮ ಅಂಕೆ_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  ವರ್ಣ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  ಪೂರ್ಣಾಙ್ಕ buffer_size;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ಯದಿ (trans == 'N') i = 0;
  ಯದಿ (trans == 'T') i = 1;
  ಯದಿ (trans == 'R') i = 0;
  ಯದಿ (trans == 'C') i = 1;

  ಯದಿ (incy == 0)	info = 11;
  ಯದಿ (incx == 0)	info = 8;
  ಯದಿ (lda < MAX(1, m))	info = 6;
  ಯದಿ (n < 0)		info = 3;
  ಯದಿ (m < 0)		info = 2;
  ಯದಿ (i < 0)          info = 1;

  trans = i;

  ಯದಿ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  ಪೂರ್ಣಾಙ್ಕ trans, buffer_size;
  blasint info, t;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 0;
    ಯದಿ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ಯದಿ (incy == 0)	  info = 11;
    ಯದಿ (incx == 0)	  info = 8;
    ಯದಿ (lda < MAX(1, m))  info = 6;
    ಯದಿ (n < 0)		  info = 3;
    ಯದಿ (m < 0)		  info = 2;
    ಯದಿ (trans < 0)        info = 1;

  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 1;
    ಯದಿ (TransA == CblasTrans)       trans = 0;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 1;
    ಯದಿ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    ಯದಿ (incy == 0)	  info = 11;
    ಯದಿ (incx == 0)	  info = 8;
    ಯದಿ (lda < MAX(1, m))  info = 6;
    ಯದಿ (n < 0)		  info = 3;
    ಯದಿ (m < 0)		  info = 2;
    ಯದಿ (trans < 0)        info = 1;

  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif
  ಯದಿ ((m==0) || (n==0)) ಪ್ರತ್ಯಾಗಮ;

  lenx = n;
  leny = m;
  ಯದಿ (trans) lenx = m;
  ಯದಿ (trans) leny = n;

  ಯದಿ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ (alpha == ZERO) ಪ್ರತ್ಯಾಗಮ;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (lenx - 1) * incx;
  ಯದಿ (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / ಪರಿಮಾಣ(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / ಪರಿಮಾಣ(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ರಾಶಿ_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ಯದಿ (nthreads == 1) {
#endif

    (gemv[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ಅನ್ಯಥಾ {

    (gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  ರಾಶಿ_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;

}
