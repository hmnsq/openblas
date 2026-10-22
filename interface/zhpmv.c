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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <ctype.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHPMV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#ଅନ୍ୟଥା
#define ERROR_NAME "CHPMV "
#endif

ସ୍ଥିର  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, ଶୂନ୍ୟ *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#ନଚେତଯଦି defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#ଅନ୍ୟଥା
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
ସ୍ଥିର  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#ନଚେତଯଦି defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#ଅନ୍ୟଥା
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ଅକ୍ଷର uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  FLOAT *buffer;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ଯଦି (uplo_arg  == 'U') uplo  = 0;
  ଯଦି (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ଯଦି (incy == 0)          info =  9;
  ଯଦି (incx == 0)          info =  6;
  ଯଦି (n < 0)              info =  2;
  ଯଦି (uplo  < 0)          info =  1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order,
	   ଗଣନା CBLAS_UPLO Uplo,
	   blasint n,
	   ଶୂନ୍ୟ *VALPHA,
	   ଶୂନ୍ୟ  *va,
	   ଶୂନ୍ୟ  *vx, blasint incx,
	   ଶୂନ୍ୟ *VBETA,
	   ଶୂନ୍ୟ  *vy, blasint incy){

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
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା uplo;
  blasint info;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ଯଦି (order == CblasColMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 0;
    ଯଦି (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ଯଦି (incy == 0)          info =  9;
    ଯଦି (incx == 0)          info =  6;
    ଯଦି (n < 0)              info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (Uplo == CblasUpper)         uplo  = 3;
    ଯଦି (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ଯଦି (incy == 0)          info =  9;
    ଯଦି (incx == 0)          info =  6;
    ଯଦି (n < 0)              info =  2;
    ଯଦି (uplo  < 0)          info =  1;
  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif

  ଯଦି (n == 0) ଫେରନ୍ତୁ;

  ଯଦି ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ଖାଲି, 0, ଖାଲି, 0);

  ଯଦି ((alpha_r == ZERO) && (alpha_i == ZERO)) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0 ) x -= (n - 1) * incx * 2;
  ଯଦି (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ଯଦି (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ifdef SMP
  nthreads = ଅଙ୍କ_cpu_avail(2);

  ଯଦି (n < 50) nthreads = 1;
  ଯଦି (nthreads > 2 && n < 500) nthreads = 2;

  ଯଦି (nthreads == 1) {
#endif

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } ଅନ୍ୟଥା {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
