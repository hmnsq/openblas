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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#അന്തര്ഭാവ "l1param.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#അല്ലെങ്കിലെങ്കിൽ defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#അല്ലെങ്കിൽ
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 ബാഹ്യ അക്ഷരം* gotoblas_corename(ശൂന്യം);
#endif

#എങ്കിൽ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, പൂർണ്ണസംഖ്യ ncpu) {
  #ifdef DOUBLE
      തിരികെ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #അല്ലെങ്കിൽ
      തിരികെ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#എങ്കിൽ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, പൂർണ്ണസംഖ്യ ncpu) {
  തിരികെ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

സ്ഥിതം വരിയിൽ പൂർണ്ണസംഖ്യ get_gemv_optimal_nthreads(BLASLONG MN) {
  പൂർണ്ണസംഖ്യ ncpu = അക്കം_cpu_avail(3);
#എങ്കിൽ defined(_WIN64) && defined(_M_ARM64)
  എങ്കിൽ (MN > 100000000L)
    തിരികെ അക്കം_cpu_avail(4);
  തിരികെ 1;
#endif
#എങ്കിൽ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  തിരികെ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#അല്ലെങ്കിലെങ്കിൽ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  തിരികെ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#അല്ലെങ്കിലെങ്കിൽ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  എങ്കിൽ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    തിരികെ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  എങ്കിൽ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    തിരികെ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  എങ്കിൽ ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    തിരികെ 1;
  അല്ലെങ്കിൽ
    തിരികെ അക്കം_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  അക്ഷരം trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  പൂർണ്ണസംഖ്യ buffer_size;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  എങ്കിൽ (trans == 'N') i = 0;
  എങ്കിൽ (trans == 'T') i = 1;
  എങ്കിൽ (trans == 'R') i = 0;
  എങ്കിൽ (trans == 'C') i = 1;

  എങ്കിൽ (incy == 0)	info = 11;
  എങ്കിൽ (incx == 0)	info = 8;
  എങ്കിൽ (lda < MAX(1, m))	info = 6;
  എങ്കിൽ (n < 0)		info = 3;
  എങ്കിൽ (m < 0)		info = 2;
  എങ്കിൽ (i < 0)          info = 1;

  trans = i;

  എങ്കിൽ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  പൂർണ്ണസംഖ്യ trans, buffer_size;
  blasint info, t;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 0;
    എങ്കിൽ (TransA == CblasTrans)       trans = 1;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    എങ്കിൽ (incy == 0)	  info = 11;
    എങ്കിൽ (incx == 0)	  info = 8;
    എങ്കിൽ (lda < MAX(1, m))  info = 6;
    എങ്കിൽ (n < 0)		  info = 3;
    എങ്കിൽ (m < 0)		  info = 2;
    എങ്കിൽ (trans < 0)        info = 1;

  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 1;
    എങ്കിൽ (TransA == CblasTrans)       trans = 0;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 1;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    എങ്കിൽ (incy == 0)	  info = 11;
    എങ്കിൽ (incx == 0)	  info = 8;
    എങ്കിൽ (lda < MAX(1, m))  info = 6;
    എങ്കിൽ (n < 0)		  info = 3;
    എങ്കിൽ (m < 0)		  info = 2;
    എങ്കിൽ (trans < 0)        info = 1;

  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif
  എങ്കിൽ ((m==0) || (n==0)) തിരികെ;

  lenx = n;
  leny = m;
  എങ്കിൽ (trans) lenx = m;
  എങ്കിൽ (trans) leny = n;

  എങ്കിൽ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ശൂന്യ, 0, ശൂന്യ, 0);

  എങ്കിൽ (alpha == ZERO) തിരികെ;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0) x -= (lenx - 1) * incx;
  എങ്കിൽ (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / വലുപ്പം(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / വലുപ്പം(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  അട്ടി_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  എങ്കിൽ (nthreads == 1) {
#endif

    (gemv[(പൂർണ്ണസംഖ്യ)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {

    (gemv_thread[(പൂർണ്ണസംഖ്യ)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  അട്ടി_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  തിരികെ;

}
