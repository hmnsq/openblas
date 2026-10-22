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
#define ERROR_NAME "QSYR2 "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DSYR2 "
#अञ्ञथा
#define ERROR_NAME "SSYR2 "
#endif

स्थिर पूर्णाङ्क (*syr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qsyr2_U, qsyr2_L,
#नोचेद्यदि defined(DOUBLE)
  dsyr2_U, dsyr2_L,
#अञ्ञथा
  ssyr2_U, ssyr2_L,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*syr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  qsyr2_thread_U, qsyr2_thread_L,
#नोचेद्यदि defined(DOUBLE)
  dsyr2_thread_U, dsyr2_thread_L,
#अञ्ञथा
  ssyr2_thread_U, ssyr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  सचे (uplo_arg  == 'U') uplo  = 0;
  सचे (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  सचे (lda  < MAX(1, n))   info =  9;
  सचे (incy == 0)          info =  7;
  सचे (incx == 0)          info =  5;
  सचे (n < 0)              info =  2;
  सचे (uplo  < 0)          info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy, FLOAT *a, blasint lda) {

  FLOAT *buffer;
  पूर्णाङ्क uplo;
  blasint info;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  सचे (order == CblasColMajor) {

    सचे (Uplo == CblasUpper) uplo  = 0;
    सचे (Uplo == CblasLower) uplo  = 1;

    info = -1;

    सचे (lda  < MAX(1, n))   info =  9;
    सचे (incy == 0)          info =  7;
    सचे (incx == 0)          info =  5;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (order == CblasRowMajor) {

    सचे (Uplo == CblasUpper) uplo  = 1;
    सचे (Uplo == CblasLower) uplo  = 0;

    info = -1;

    सचे (lda  < MAX(1, n))   info =  9;
    सचे (incy == 0)          info =  7;
    सचे (incx == 0)          info =  5;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  सचे (n == 0) पटिनिवत्त;

  सचे (alpha == ZERO) पटिनिवत्त;

  IDEBUG_START;

  सचे (incx == 1 && incy == 1 && n < 100) {
    blasint i;
    सचे (!uplo) {
      अत्थाय (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, रिक्त, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, रिक्त, 0);
        a += lda;
      }
    } अञ्ञथा {
      अत्थाय (i = 0; i < n; i++){
        AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, रिक्त, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, रिक्त, 0);
        a += 1 + lda;
      }
    }
    पटिनिवत्त;
  }

	  
  FUNCTION_PROFILE_START();

  सचे (incx < 0 ) x -= (n - 1) * incx;
  सचे (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  सचे (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#ifdef SMP
  nthreads = अङ्कः_cpu_avail(2);

  सचे (n < 50) nthreads = 1;
  सचे (nthreads > 2 && n < 500) nthreads = 2;

  सचे (nthreads == 1) {
#endif

    (syr2[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer);

#ifdef SMP
  } अञ्ञथा {

    (syr2_thread[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  पटिनिवत्त;
}
