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
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMV "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMV "
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMV "
#اگر_ختم

#اگر_تعریف_شدہ SMP
ساکن صحیح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#اگر_تعریف_شدہ XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#وگرنہ_اگر defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ورنہ
  sgemv_thread_n, sgemv_thread_t,
#اگر_ختم
};
#اگر_ختم

#اگر_تعریف_شدہ SMP
#اگر_تعریف_شدہ DYNAMIC_ARCH
 خارجی حرف* gotoblas_corename(خالی);
#اگر_ختم

#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ساکن اندرونی صحیح get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, صحیح ncpu) {
  #اگر_تعریف_شدہ DOUBLE
      واپس (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ورنہ
      واپس
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #اگر_ختم
}
#اگر_ختم

#اگر defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ساکن اندرونی صحیح get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, صحیح ncpu) {
  واپس
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#اگر_ختم

ساکن اندرونی صحیح get_gemv_optimal_nthreads(BLASLONG MN) {
  صحیح ncpu = عدد_cpu_avail(3);
#اگر defined(_WIN64) && defined(_M_ARM64)
  اگر (MN > 100000000L)
    واپس عدد_cpu_avail(4);
  واپس 1;
#اگر_ختم
#اگر defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  واپس get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#وگرنہ_اگر defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  واپس get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#وگرنہ_اگر defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  اگر (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    واپس get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  اگر (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    واپس get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#اگر_ختم

  اگر ( MN < 115200L * GEMM_MULTITHREAD_حد )
    واپس 1;
  ورنہ
    واپس عدد_cpu_avail(2);
}
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *TRANS, blasint *M, blasint *N,
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
  FLOAT *بفر;
  صحیح بفر_سائز;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

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
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب,
	   شمار CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *بفر;
  blasint lenx, leny;
  صحیح trans, بفر_سائز;
  blasint info, t;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  صحیح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
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

  اگر (ترتیب == Cblas_صف_Major) {
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
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_ختم
  اگر ((m==0) || (n==0)) واپس;

  lenx = n;
  leny = m;
  اگر (trans) lenx = m;
  اگر (trans) leny = n;

  اگر (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), تہی, 0, تہی, 0);

  اگر (alpha == ZERO) واپس;
	
  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0) x -= (lenx - 1) * incx;
  اگر (incy < 0) y -= (leny - 1) * incy;

  بفر_سائز = m + n + 128 / حجم(FLOAT);
#اگر_تعریف_شدہ WINDOWS_ABI
  بفر_سائز += 160 / حجم(FLOAT) ;
#اگر_ختم
  // for alignment
  بفر_سائز = (بفر_سائز + 3) & ~3;
  ڈھیر_ALLOC(بفر_سائز, FLOAT, بفر);

#اگر_تعریف_شدہ SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  اگر (nthreads == 1) {
#اگر_ختم

    (gemv[(صحیح)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, بفر);

#اگر_تعریف_شدہ SMP
  } ورنہ {

    (gemv_thread[(صحیح)trans])(m, n, alpha, a, lda, x, incx, y, incy, بفر, nthreads);

  }
#اگر_ختم

  ڈھیر_FREE(بفر);
  FUNCTION_PROFILE_انتہا(1, m * n + m + n,  2 * m * n);

  IDEBUG_انتہا;

  واپس;

}
