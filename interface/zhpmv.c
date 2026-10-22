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
#define ERROR_NAME "XHPMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#अञ्ञथा
#define ERROR_NAME "CHPMV "
#endif

स्थिर  पूर्णाङ्क (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#नोचेद्यदि defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#अञ्ञथा
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
स्थिर  पूर्णाङ्क (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#नोचेद्यदि defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#अञ्ञथा
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
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

  सचे (incy == 0)          info =  9;
  सचे (incx == 0)          info =  6;
  सचे (n < 0)              info =  2;
  सचे (uplo  < 0)          info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_UPLO Uplo,
	   blasint n,
	   शून्य *VALPHA,
	   शून्य  *va,
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
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  सचे (order == CblasColMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 0;
    सचे (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    सचे (incy == 0)          info =  9;
    सचे (incx == 0)          info =  6;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (order == CblasRowMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 3;
    सचे (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    सचे (incy == 0)          info =  9;
    सचे (incx == 0)          info =  6;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  सचे (n == 0) पटिनिवत्त;

  सचे ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  सचे ((alpha_r == ZERO) && (alpha_i == ZERO)) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0 ) x -= (n - 1) * incx * 2;
  सचे (incy < 0 ) y -= (n - 1) * incy * 2;

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

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } अञ्ञथा {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  पटिनिवत्त;
}
