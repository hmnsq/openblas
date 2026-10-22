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
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DGEMV "
#अन्यथा
#define ERROR_NAME "SGEMV "
#endif

#ifdef SMP
स्थिर पूर्णाङ्क (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#नोचेद्यदि defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#अन्यथा
  sgemv_thread_n, sgemv_thread_t,
#endif
};
#endif

#ifdef SMP
#ifdef DYNAMIC_ARCH
 बाह्य वर्ण* gotoblas_corename(शून्य);
#endif

#यदि defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
स्थिर अन्तःस्थ पूर्णाङ्क get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, पूर्णाङ्क ncpu) {
  #ifdef DOUBLE
      प्रत्यागम (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #अन्यथा
      प्रत्यागम
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #endif
}
#endif

#यदि defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
स्थिर अन्तःस्थ पूर्णाङ्क get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, पूर्णाङ्क ncpu) {
  प्रत्यागम
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#endif

स्थिर अन्तःस्थ पूर्णाङ्क get_gemv_optimal_nthreads(BLASLONG MN) {
  पूर्णाङ्क ncpu = अंक_cpu_avail(3);
#यदि defined(_WIN64) && defined(_M_ARM64)
  यदि (MN > 100000000L)
    प्रत्यागम अंक_cpu_avail(4);
  प्रत्यागम 1;
#endif
#यदि defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  प्रत्यागम get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#नोचेद्यदि defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  प्रत्यागम get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#नोचेद्यदि defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  यदि (म_भेद(gotoblas_corename(), "neoversev1") == 0) {
    प्रत्यागम get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  यदि (म_भेद(gotoblas_corename(), "neoversev2") == 0) {
    प्रत्यागम get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#endif

  यदि ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    प्रत्यागम 1;
  अन्यथा
    प्रत्यागम अंक_cpu_avail(2);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  वर्ण trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *buffer;
  पूर्णाङ्क buffer_size;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  यदि (trans == 'N') i = 0;
  यदि (trans == 'T') i = 1;
  यदि (trans == 'R') i = 0;
  यदि (trans == 'C') i = 1;

  यदि (incy == 0)	info = 11;
  यदि (incx == 0)	info = 8;
  यदि (lda < MAX(1, m))	info = 6;
  यदि (n < 0)		info = 3;
  यदि (m < 0)		info = 2;
  यदि (i < 0)          info = 1;

  trans = i;

  यदि (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
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
  पूर्णाङ्क trans, buffer_size;
  blasint info, t;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  यदि (order == CblasColMajor) {
    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 0;
    यदि (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    यदि (incy == 0)	  info = 11;
    यदि (incx == 0)	  info = 8;
    यदि (lda < MAX(1, m))  info = 6;
    यदि (n < 0)		  info = 3;
    यदि (m < 0)		  info = 2;
    यदि (trans < 0)        info = 1;

  }

  यदि (order == CblasRowMajor) {
    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 1;
    यदि (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    यदि (incy == 0)	  info = 11;
    यदि (incx == 0)	  info = 8;
    यदि (lda < MAX(1, m))  info = 6;
    यदि (n < 0)		  info = 3;
    यदि (m < 0)		  info = 2;
    यदि (trans < 0)        info = 1;

  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#endif
  यदि ((m==0) || (n==0)) प्रत्यागम;

  lenx = n;
  leny = m;
  यदि (trans) lenx = m;
  यदि (trans) leny = n;

  यदि (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  यदि (alpha == ZERO) प्रत्यागम;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (lenx - 1) * incx;
  यदि (incy < 0) y -= (leny - 1) * incy;

  buffer_size = m + n + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / परिमाण(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ढेर_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  यदि (nthreads == 1) {
#endif

    (gemv[(पूर्णाङ्क)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अन्यथा {

    (gemv_thread[(पूर्णाङ्क)trans])(m, n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  ढेर_FREE(buffer);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  प्रत्यागम;

}
