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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#અન્તર્ભાવ "l1param.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#નહીંતર
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
સ્થિર પૂર્ણાંક (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#અથવાજો defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#નહીંતર
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 બાહ્ય અક્ષર* gotoblas_corename(શૂન્ય);
#endif

#જો defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
સ્થિર અંતઃરેખ પૂર્ણાંક get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, પૂર્ણાંક ncpu) {
  #ifdef DOUBLE
      પાછા (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #નહીંતર
      પાછા
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#જો defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
સ્થિર અંતઃરેખ પૂર્ણાંક get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, પૂર્ણાંક ncpu) {
  પાછા
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

સ્થિર અંતઃરેખ પૂર્ણાંક get_gemv_optimal_nthreads(BLASLONG MN) {
  પૂર્ણાંક ncpu = અંક_cpu_avail(3);
#જો defined(_WIN64) && defined(_M_ARM64)
  જો (MN > 100000000L)
    પાછા અંક_cpu_avail(4);
  પાછા 1;
#endif
#જો defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  પાછા get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#અથવાજો defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  પાછા get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#અથવાજો defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  જો (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    પાછા get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  જો (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    પાછા get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  જો ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    પાછા 1;
  નહીંતર
    પાછા અંક_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  અક્ષર trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  પૂર્ણાંક buffer_size;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  જો (trans == 'N') i = 0;
  જો (trans == 'T') i = 1;
  જો (trans == 'R') i = 0;
  જો (trans == 'C') i = 1;

  જો (incy == 0)	info = 11;
  જો (incx == 0)	info = 8;
  જો (lda < MAX(1, m))	info = 6;
  જો (n < 0)		info = 3;
  જો (m < 0)		info = 2;
  જો (i < 0)          info = 1;

  trans = i;

  જો (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   ગણના CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  પૂર્ણાંક trans, buffer_size;
  blasint info, t;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (TransA == CblasNoTrans)     trans = 0;
    જો (TransA == CblasTrans)       trans = 1;
    જો (TransA == CblasConjNoTrans) trans = 0;
    જો (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    જો (incy == 0)	  info = 11;
    જો (incx == 0)	  info = 8;
    જો (lda < MAX(1, m))  info = 6;
    જો (n < 0)		  info = 3;
    જો (m < 0)		  info = 2;
    જો (trans < 0)        info = 1;

  }

  જો (order == CblasRowMajor) {
    જો (TransA == CblasNoTrans)     trans = 1;
    જો (TransA == CblasTrans)       trans = 0;
    જો (TransA == CblasConjNoTrans) trans = 1;
    જો (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    જો (incy == 0)	  info = 11;
    જો (incx == 0)	  info = 8;
    જો (lda < MAX(1, m))  info = 6;
    જો (n < 0)		  info = 3;
    જો (m < 0)		  info = 2;
    જો (trans < 0)        info = 1;

  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif
  જો ((m==0) || (n==0)) પાછા;

  lenx = n;
  leny = m;
  જો (trans) lenx = m;
  જો (trans) leny = n;

  જો (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ખાલી, 0, ખાલી, 0);

  જો (alpha == ZERO) પાછા;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (lenx - 1) * incx;
  જો (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / કદ(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / કદ(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ઢગલો_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  જો (nthreads == 1) {
#endif

    (gemv[(પૂર્ણાંક)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } નહીંતર {

    (gemv_thread[(પૂર્ણાંક)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  ઢગલો_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  પાછા;

}
