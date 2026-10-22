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
#නැතහොත්
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
static පූර්ණ_සංඛ්‍යා (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, පූර්ණ_සංඛ්‍යා) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#elif defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#නැතහොත්
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 extern var_na* gotoblas_corename(void);
#endif

#නම් defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
static anta_hstha පූර්ණ_සංඛ්‍යා get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, පූර්ණ_සංඛ්‍යා ncpu) {
  #ifdef DOUBLE
      ආපසු (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #නැතහොත්
      ආපසු
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#නම් defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
static anta_hstha පූර්ණ_සංඛ්‍යා get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, පූර්ණ_සංඛ්‍යා ncpu) {
  ආපසු
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

static anta_hstha පූර්ණ_සංඛ්‍යා get_gemv_optimal_nthreads(BLASLONG MN) {
  පූර්ණ_සංඛ්‍යා ncpu = num_cpu_avail(3);
#නම් defined(_WIN64) && defined(_M_ARM64)
  නම් (MN > 100000000L)
    ආපසු num_cpu_avail(4);
  ආපසු 1;
#endif
#නම් defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ආපසු get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#elif defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ආපසු get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#elif defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  නම් (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ආපසු get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  නම් (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ආපසු get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  නම් ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    ආපසු 1;
  නැතහොත්
    ආපසු num_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  පූර්ණ_සංඛ්‍යා buffer_size;
#ifdef SMP
  පූර්ණ_සංඛ්‍යා nthreads;
#endif

  පූර්ණ_සංඛ්‍යා (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  නම් (trans == 'N') i = 0;
  නම් (trans == 'T') i = 1;
  නම් (trans == 'R') i = 0;
  නම් (trans == 'C') i = 1;

  නම් (incy == 0)	info = 11;
  නම් (incx == 0)	info = 8;
  නම් (lda < MAX(1, m))	info = 6;
  නම් (n < 0)		info = 3;
  නම් (m < 0)		info = 2;
  නම් (i < 0)          info = 1;

  trans = i;

  නම් (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ආපසු;
  }

#නැතහොත්

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  පූර්ණ_සංඛ්‍යා trans, buffer_size;
  blasint info, t;
#ifdef SMP
  පූර්ණ_සංඛ්‍යා nthreads;
#endif

  පූර්ණ_සංඛ්‍යා (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  නම් (order == CblasColMajor) {
    නම් (TransA == CblasNoTrans)     trans = 0;
    නම් (TransA == CblasTrans)       trans = 1;
    නම් (TransA == CblasConjNoTrans) trans = 0;
    නම් (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    නම් (incy == 0)	  info = 11;
    නම් (incx == 0)	  info = 8;
    නම් (lda < MAX(1, m))  info = 6;
    නම් (n < 0)		  info = 3;
    නම් (m < 0)		  info = 2;
    නම් (trans < 0)        info = 1;

  }

  නම් (order == CblasRowMajor) {
    නම් (TransA == CblasNoTrans)     trans = 1;
    නම් (TransA == CblasTrans)       trans = 0;
    නම් (TransA == CblasConjNoTrans) trans = 1;
    නම් (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    නම් (incy == 0)	  info = 11;
    නම් (incx == 0)	  info = 8;
    නම් (lda < MAX(1, m))  info = 6;
    නම් (n < 0)		  info = 3;
    නම් (m < 0)		  info = 2;
    නම් (trans < 0)        info = 1;

  }

  නම් (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ආපසු;
  }

#endif
  නම් ((m==0) || (n==0)) ආපසු;

  lenx = n;
  leny = m;
  නම් (trans) lenx = m;
  නම් (trans) leny = n;

  නම් (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  නම් (alpha == ZERO) ආපසු;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  නම් (incx < 0) x -= (lenx - 1) * incx;
  නම් (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / parimaa_na(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / parimaa_na(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  නම් (nthreads == 1) {
#endif

    (gemv[(පූර්ණ_සංඛ්‍යා)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } නැතහොත් {

    (gemv_thread[(පූර්ණ_සංඛ්‍යා)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  STACK_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ආපසු;

}
