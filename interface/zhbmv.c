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
#define ERROR_NAME "XHBMV "
#अथवायदि defined(DOUBLE)
#define ERROR_NAME "ZHBMV "
#अन्यथा
#define ERROR_NAME "CHBMV "
#endif

स्थिर  पूर्णाङ्क (*hbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  xhbmv_U, xhbmv_L, xhbmv_V, xhbmv_M,
#अथवायदि defined(DOUBLE)
  zhbmv_U, zhbmv_L, zhbmv_V, zhbmv_M,
#अन्यथा
  chbmv_U, chbmv_L, chbmv_V, chbmv_M,
#endif
};

#ifdef SMPBUG
स्थिर  पूर्णाङ्क (*hbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xhbmv_thread_U, xhbmv_thread_L, xhbmv_thread_V, xhbmv_thread_M,
#अथवायदि defined(DOUBLE)
  zhbmv_thread_U, zhbmv_thread_L, zhbmv_thread_V, zhbmv_thread_M,
#अन्यथा
  chbmv_thread_U, chbmv_thread_L, chbmv_thread_V, chbmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  वर्ण uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *buffer;
#ifdef SMPBUG
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;
  यदि (uplo_arg  == 'V') uplo  = 2;
  यदि (uplo_arg  == 'M') uplo  = 3;

  info = 0;

  यदि (incy == 0)          info = 11;
  यदि (incx == 0)          info =  8;
  यदि (lda  < k + 1)       info =  6;
  यदि (k < 0)              info =  3;
  यदि (n < 0)              info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   शून्य *VALPHA,
	   शून्य  *va, blasint lda,
	   शून्य  *vx, blasint incx,
	   शून्य *VBETA,
	   शून्य  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *buffer;
  पूर्णाङ्क uplo;
  blasint info;
#ifdef SMPBUG
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  यदि (order == CblasColMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 0;
    यदि (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    यदि (incy == 0)          info = 11;
    यदि (incx == 0)          info =  8;
    यदि (lda  < k + 1)       info =  6;
    यदि (k < 0)              info =  3;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (order == CblasRowMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 3;
    यदि (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    यदि (incy == 0)          info = 11;
    यदि (incx == 0)          info =  8;
    यदि (lda  < k + 1)       info =  6;
    यदि (k < 0)              info =  3;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#endif

  यदि (n == 0) फर्काऊ;

  यदि ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  यदि ((alpha_r == ZERO) && (alpha_i == ZERO)) फर्काऊ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  यदि (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  यदि (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#ifdef SMPBUG
  nthreads = अङ्क_cpu_avail(2);

  यदि (nthreads == 1) {
#endif

    (hbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMPBUG
  } अन्यथा {

    (hbmv_thread[uplo])(n, k, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n,  n * k);

  IDEBUG_END;

  फर्काऊ;
}
