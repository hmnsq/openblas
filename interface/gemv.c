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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#अन्तर्भाव "l1param.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGEMV "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#अन्यथा
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
स्थिर पूर्णांक (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#किंवाजर defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#अन्यथा
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 बाह्य अक्षर* gotoblas_corename(शून्य);
#endif

#जर defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
स्थिर ओळीत पूर्णांक get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, पूर्णांक ncpu) {
  #ifdef DOUBLE
      परत (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #अन्यथा
      परत
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#जर defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
स्थिर ओळीत पूर्णांक get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, पूर्णांक ncpu) {
  परत
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

स्थिर ओळीत पूर्णांक get_gemv_optimal_nthreads(BLASLONG MN) {
  पूर्णांक ncpu = अंक_cpu_avail(3);
#जर defined(_WIN64) && defined(_M_ARM64)
  जर (MN > 100000000L)
    परत अंक_cpu_avail(4);
  परत 1;
#endif
#जर defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  परत get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#किंवाजर defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  परत get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#किंवाजर defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  जर (म_भेद(gotoblas_corename(), "neoversev1") == 0) {
    परत get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  जर (म_भेद(gotoblas_corename(), "neoversev2") == 0) {
    परत get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  जर ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    परत 1;
  अन्यथा
    परत अंक_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  पूर्णांक buffer_size;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  जर (trans == 'N') i = 0;
  जर (trans == 'T') i = 1;
  जर (trans == 'R') i = 0;
  जर (trans == 'C') i = 1;

  जर (incy == 0)	info = 11;
  जर (incx == 0)	info = 8;
  जर (lda < MAX(1, m))	info = 6;
  जर (n < 0)		info = 3;
  जर (m < 0)		info = 2;
  जर (i < 0)          info = 1;

  trans = i;

  जर (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny;
  पूर्णांक trans, buffer_size;
  blasint info, t;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  जर (order == CblasColMajor) {
    जर (TransA == CblasNoTrans)     trans = 0;
    जर (TransA == CblasTrans)       trans = 1;
    जर (TransA == CblasConjNoTrans) trans = 0;
    जर (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    जर (incy == 0)	  info = 11;
    जर (incx == 0)	  info = 8;
    जर (lda < MAX(1, m))  info = 6;
    जर (n < 0)		  info = 3;
    जर (m < 0)		  info = 2;
    जर (trans < 0)        info = 1;

  }

  जर (order == CblasRowMajor) {
    जर (TransA == CblasNoTrans)     trans = 1;
    जर (TransA == CblasTrans)       trans = 0;
    जर (TransA == CblasConjNoTrans) trans = 1;
    जर (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    जर (incy == 0)	  info = 11;
    जर (incx == 0)	  info = 8;
    जर (lda < MAX(1, m))  info = 6;
    जर (n < 0)		  info = 3;
    जर (m < 0)		  info = 2;
    जर (trans < 0)        info = 1;

  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#endif
  जर ((m==0) || (n==0)) परत;

  lenx = n;
  leny = m;
  जर (trans) lenx = m;
  जर (trans) leny = n;

  जर (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  जर (alpha == ZERO) परत;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जर (incx < 0) x -= (lenx - 1) * incx;
  जर (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / आकार(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / आकार(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  चळत_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  जर (nthreads == 1) {
#endif

    (gemv[(पूर्णांक)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अन्यथा {

    (gemv_thread[(पूर्णांक)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  चळत_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  परत;

}
