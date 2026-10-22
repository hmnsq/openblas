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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QGEMV "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DGEMV "
#وإلا
#وصف_ڪريو ERROR_NAME "SGEMV "
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
ساكن صحيح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#وإلاإذا defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#وإلا
  sgemv_thread_n, sgemv_thread_t,
#جيڪڏهن_پورو
};
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_ٿيل DYNAMIC_ARCH
 خارجي حرف* gotoblas_corename(فراغ);
#جيڪڏهن_پورو

#إذا defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ساكن مضمن صحيح get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, صحيح ncpu) {
  #جيڪڏهن_وصف_ٿيل DOUBLE
      أرجع (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #وإلا
      أرجع
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #جيڪڏهن_پورو
}
#جيڪڏهن_پورو

#إذا defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ساكن مضمن صحيح get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, صحيح ncpu) {
  أرجع
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#جيڪڏهن_پورو

ساكن مضمن صحيح get_gemv_optimal_nthreads(BLASLONG MN) {
  صحيح ncpu = انگ_cpu_avail(3);
#إذا defined(_WIN64) && defined(_M_ARM64)
  إذا (MN > 100000000L)
    أرجع انگ_cpu_avail(4);
  أرجع 1;
#جيڪڏهن_پورو
#إذا defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  أرجع get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#وإلاإذا defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  أرجع get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#وإلاإذا defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  إذا (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    أرجع get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  إذا (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    أرجع get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#جيڪڏهن_پورو

  إذا ( MN < 115200L * GEMM_MULTITHREAD_حد )
    أرجع 1;
  وإلا
    أرجع انگ_cpu_avail(2);
}
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

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
  FLOAT *بفر;
  صحيح بفر_سائيز;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

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

  إذا (trans == 'N') i = 0;
  إذا (trans == 'T') i = 1;
  إذا (trans == 'R') i = 0;
  إذا (trans == 'C') i = 1;

  إذا (incy == 0)	info = 11;
  إذا (incx == 0)	info = 8;
  إذا (lda < MAX(1, m))	info = 6;
  إذا (n < 0)		info = 3;
  إذا (m < 0)		info = 2;
  إذا (i < 0)          info = 1;

  trans = i;

  إذا (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *بفر;
  blasint lenx, leny;
  صحيح trans, بفر_سائيز;
  blasint info, t;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  صحيح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 0;
    إذا (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    إذا (incy == 0)	  info = 11;
    إذا (incx == 0)	  info = 8;
    إذا (lda < MAX(1, m))  info = 6;
    إذا (n < 0)		  info = 3;
    إذا (m < 0)		  info = 2;
    إذا (trans < 0)        info = 1;

  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 1;
    إذا (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    إذا (incy == 0)	  info = 11;
    إذا (incx == 0)	  info = 8;
    إذا (lda < MAX(1, m))  info = 6;
    إذا (n < 0)		  info = 3;
    إذا (m < 0)		  info = 2;
    إذا (trans < 0)        info = 1;

  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_پورو
  إذا ((m==0) || (n==0)) أرجع;

  lenx = n;
  leny = m;
  إذا (trans) lenx = m;
  إذا (trans) leny = n;

  إذا (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  إذا (alpha == ZERO) أرجع;
	
  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (lenx - 1) * incx;
  إذا (incy < 0) y -= (leny - 1) * incy;

  بفر_سائيز = m + n + 128 / حجم(FLOAT);
#جيڪڏهن_وصف_ٿيل WINDOWS_ABI
  بفر_سائيز += 160 / حجم(FLOAT) ;
#جيڪڏهن_پورو
  // for alignment
  بفر_سائيز = (بفر_سائيز + 3) & ~3;
  ڍير_ALLOC(بفر_سائيز, FLOAT, بفر);

#جيڪڏهن_وصف_ٿيل SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

    (gemv[(صحيح)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, بفر);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

    (gemv_thread[(صحيح)trans])(m, n, alpha, a, lda, x, incx, y, incy, بفر, nthreads);

  }
#جيڪڏهن_پورو

  ڍير_FREE(بفر);
  FUNCTION_PROFILE_پڇاڙي(1, m * n + m + n,  2 * m * n);

  IDEBUG_پڇاڙي;

  أرجع;

}
