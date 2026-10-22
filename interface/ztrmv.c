/*********************************************************************/
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
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XTRMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZTRMV "
#अण्णहा
#define ERROR_NAME "CTRMV "
#endif

स्थिर पूर्णाङ्क (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_RUN, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_RUN, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#अण्णहा
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_RUN, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_RUN, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_RUN, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#अण्णहा
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_RUN, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans, buffer_size;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  जइ (trans_arg == 'N') trans = 0;
  जइ (trans_arg == 'T') trans = 1;
  जइ (trans_arg == 'R') trans = 2;
  जइ (trans_arg == 'C') trans = 3;

  जइ (diag_arg  == 'U') unit  = 0;
  जइ (diag_arg  == 'N') unit  = 1;

  जइ (uplo_arg  == 'U') uplo  = 0;
  जइ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  जइ (incx == 0)          info =  8;
  जइ (lda  < MAX(1, n))   info =  6;
  जइ (n < 0)              info =  4;
  जइ (unit  < 0)          info =  3;
  जइ (trans < 0)          info =  2;
  जइ (uplo  < 0)          info =  1;

  जइ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अण्णहा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, शून्य  *va, blasint lda, शून्य  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  पूर्णाङ्क trans, uplo, unit, buffer_size;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  जइ (order == CblasColMajor) {
    जइ (Uplo == CblasUpper)         uplo  = 0;
    जइ (Uplo == CblasLower)         uplo  = 1;

    जइ (TransA == CblasNoTrans)     trans = 0;
    जइ (TransA == CblasTrans)       trans = 1;
    जइ (TransA == CblasConjNoTrans) trans = 2;
    जइ (TransA == CblasConjTrans)   trans = 3;

    जइ (Diag == CblasUnit)          unit  = 0;
    जइ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    जइ (incx == 0)          info =  8;
    जइ (lda  < MAX(1, n))   info =  6;
    जइ (n < 0)              info =  4;
    जइ (unit  < 0)          info =  3;
    जइ (trans < 0)          info =  2;
    जइ (uplo  < 0)          info =  1;
  }

  जइ (order == CblasRowMajor) {
    जइ (Uplo == CblasUpper)         uplo  = 1;
    जइ (Uplo == CblasLower)         uplo  = 0;

    जइ (TransA == CblasNoTrans)     trans = 1;
    जइ (TransA == CblasTrans)       trans = 0;
    जइ (TransA == CblasConjNoTrans) trans = 3;
    जइ (TransA == CblasConjTrans)   trans = 2;

    जइ (Diag == CblasUnit)          unit  = 0;
    जइ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    जइ (incx == 0)          info =  8;
    जइ (lda  < MAX(1, n))   info =  6;
    जइ (n < 0)              info =  4;
    जइ (unit  < 0)          info =  3;
    जइ (trans < 0)          info =  2;
    जइ (uplo  < 0)          info =  1;
  }

  जइ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#endif

  जइ (n == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जइ (incx < 0 ) x -= (n - 1) * incx * 2;

#ifdef SMP
  // Calibrated on a Xeon E5-2630
  जइ(1L * n * n > 36L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_THRESHOLD) {
    nthreads = अङ्कः_cpu_avail(2);
    जइ(nthreads > 2 && 1L * n * n < 64L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_THRESHOLD)
      nthreads = 2;
  } अण्णहा
      nthreads = 1;

  जइ(nthreads > 1) {
    buffer_size = n > 16 ? 0 : n * 4 + 40;
  }
  अण्णहा
#endif
  {
    buffer_size = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / परिमाण(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    buffer_size += 8;
    जइ(incx != 1)
      buffer_size += n * 2;
  }
  राशिः_ALLOC(buffer_size, FLOAT, buffer);

#ifdef SMP
  जइ (nthreads == 1) {
#endif

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

#ifdef SMP
  } अण्णहा {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer, nthreads);

  }
#endif

  राशिः_FREE(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n,  n * n);

  IDEBUG_END;

  प्रत्यागम;
}
