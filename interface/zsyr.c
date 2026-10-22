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
#define ERROR_NAME "XSYR  "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYR  "
#अञ्ञथा
#define ERROR_NAME "CSYR  "
#endif

स्थिर पूर्णाङ्क (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr_U, xsyr_L,
#नोचेद्यदि defined(DOUBLE)
  zsyr_U, zsyr_L,
#अञ्ञथा
  csyr_U, csyr_L,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#नोचेद्यदि defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#अञ्ञथा
  csyr_thread_U, csyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

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

  सचे (lda  < MAX(1, n))   info =  7;
  सचे (incx == 0)          info =  5;
  सचे (n < 0)              info =  2;
  सचे (uplo  < 0)          info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }


#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, पूर्णाङ्क n, शून्य* valpha, FLOAT *x, पूर्णाङ्क incx, FLOAT *a, पूर्णाङ्क lda) {

  FLOAT *buffer;
  पूर्णाङ्क uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
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

    सचे (lda  < MAX(1, n))   info =  7;
    सचे (incx == 0)          info =  5;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;

  }

  सचे (order == CblasRowMajor) {

    सचे (Uplo == CblasUpper) uplo  = 1;
    सचे (Uplo == CblasLower) uplo  = 0;

    info = -1;

    सचे (lda  < MAX(1, n))   info =  7;
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

  सचे ((alpha_r == ZERO) && (alpha_i == ZERO)) पटिनिवत्त;

  सचे (incx == 1 && n < 50) {
    blasint i;
    सचे (!uplo) {
      अत्थाय (i = 0; i < n; i++){
        सचे ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, रिक्त, 0);
        }
        a += lda * 2;
      }
    } अञ्ञथा {
      अत्थाय (i = 0; i < n; i++){
        सचे ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, रिक्त, 0);
        }
        a += 2 + lda * 2;
      }
    }
    पटिनिवत्त;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0 ) x -= (n - 1) * incx * 2;

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

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, buffer);

#ifdef SMP
  } अञ्ञथा {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  पटिनिवत्त;
}
