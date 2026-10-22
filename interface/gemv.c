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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਅਨ੍ਤਰ੍ਭਾਵ "l1param.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGEMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGEMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGEMV "
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#ਜਾਂਜੇ defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ਨਹੀਂਤਾਂ
  sgemv_thread_n, sgemv_thread_t,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
 ਬਾਹਰੀ ਅੱਖਰ* gotoblas_corename(ਖਾਲੀ);
#ਜੇ_ਅੰਤ

#ਜੇ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, ਪੂਰਨ ncpu) {
  #ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
      ਵਾਪਸ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ਨਹੀਂਤਾਂ
      ਵਾਪਸ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #ਜੇ_ਅੰਤ
}
#ਜੇ_ਅੰਤ

#ਜੇ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, ਪੂਰਨ ncpu) {
  ਵਾਪਸ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਸਤਰਅੰਦਰ ਪੂਰਨ get_gemv_optimal_nthreads(BLASLONG MN) {
  ਪੂਰਨ ncpu = ਅੰਕ_cpu_avail(3);
#ਜੇ defined(_WIN64) && defined(_M_ARM64)
  ਜੇ (MN > 100000000L)
    ਵਾਪਸ ਅੰਕ_cpu_avail(4);
  ਵਾਪਸ 1;
#ਜੇ_ਅੰਤ
#ਜੇ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ਵਾਪਸ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#ਜਾਂਜੇ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ਵਾਪਸ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#ਜਾਂਜੇ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ਜੇ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ਵਾਪਸ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ਜੇ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ਵਾਪਸ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#ਜੇ_ਅੰਤ

  ਜੇ ( MN < 115200L * GEMM_MULTITHREAD_ਹੱਦ )
    ਵਾਪਸ 1;
  ਨਹੀਂਤਾਂ
    ਵਾਪਸ ਅੰਕ_cpu_avail(2);
}
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  ਅੱਖਰ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *ਬਫਰ;
  ਪੂਰਨ ਬਫਰ_ਮਾਪ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ਜੇ (trans == 'N') i = 0;
  ਜੇ (trans == 'T') i = 1;
  ਜੇ (trans == 'R') i = 0;
  ਜੇ (trans == 'C') i = 1;

  ਜੇ (incy == 0)	info = 11;
  ਜੇ (incx == 0)	info = 8;
  ਜੇ (lda < MAX(1, m))	info = 6;
  ਜੇ (n < 0)		info = 3;
  ਜੇ (m < 0)		info = 2;
  ਜੇ (i < 0)          info = 1;

  trans = i;

  ਜੇ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *ਬਫਰ;
  blasint lenx, leny;
  ਪੂਰਨ trans, ਬਫਰ_ਮਾਪ;
  blasint info, t;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 0;
    ਜੇ (TransA == CblasTrans)       trans = 1;
    ਜੇ (TransA == CblasConjNoTrans) trans = 0;
    ਜੇ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ਜੇ (incy == 0)	  info = 11;
    ਜੇ (incx == 0)	  info = 8;
    ਜੇ (lda < MAX(1, m))  info = 6;
    ਜੇ (n < 0)		  info = 3;
    ਜੇ (m < 0)		  info = 2;
    ਜੇ (trans < 0)        info = 1;

  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 1;
    ਜੇ (TransA == CblasTrans)       trans = 0;
    ਜੇ (TransA == CblasConjNoTrans) trans = 1;
    ਜੇ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    ਜੇ (incy == 0)	  info = 11;
    ਜੇ (incx == 0)	  info = 8;
    ਜੇ (lda < MAX(1, m))  info = 6;
    ਜੇ (n < 0)		  info = 3;
    ਜੇ (m < 0)		  info = 2;
    ਜੇ (trans < 0)        info = 1;

  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ
  ਜੇ ((m==0) || (n==0)) ਵਾਪਸ;

  lenx = n;
  leny = m;
  ਜੇ (trans) lenx = m;
  ਜੇ (trans) leny = n;

  ਜੇ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);

  ਜੇ (alpha == ZERO) ਵਾਪਸ;
	
  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (lenx - 1) * incx;
  ਜੇ (incy < 0) y -= (leny - 1) * incy;

  ਬਫਰ_ਮਾਪ = m + n + 128 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ WINDOWS_ABI
  ਬਫਰ_ਮਾਪ += 160 / ਆਕਾਰ(FLOAT) ;
#ਜੇ_ਅੰਤ
  // for alignment
  ਬਫਰ_ਮਾਪ = (ਬਫਰ_ਮਾਪ + 3) & ~3;
  ਢੇਰ_ALLOC(ਬਫਰ_ਮਾਪ, FLOAT, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

    (gemv[(ਪੂਰਨ)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    (gemv_thread[(ਪੂਰਨ)trans])(m, n, alpha, a, lda, x, incx, y, incy, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  ਢੇਰ_FREE(ਬਫਰ);
  FUNCTION_PROFILE_ਸਿਰਾ(1, m * n + m + n,  2 * m * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;

}
