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

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHPMV "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#وگرنه
#define ERROR_NAME "CHPMV "
#endif

ایستا  صحیح (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, تهی *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#وگرنهاگر defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#وگرنه
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
ایستا  صحیح (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#وگرنهاگر defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#وگرنه
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  نویسه uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  صحیح uplo;
  FLOAT *buffer;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  اگر (incy == 0)          info =  9;
  اگر (incx == 0)          info =  6;
  اگر (n < 0)              info =  2;
  اگر (uplo  < 0)          info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order,
	   شمارش CBLAS_UPLO Uplo,
	   blasint n,
	   تهی *VALPHA,
	   تهی  *va,
	   تهی  *vx, blasint incx,
	   تهی *VBETA,
	   تهی  *vy, blasint incy){

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
  صحیح uplo;
  blasint info;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  اگر (order == CblasColMajor) {
    اگر (Uplo == CblasUpper)         uplo  = 0;
    اگر (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    اگر (incy == 0)          info =  9;
    اگر (incx == 0)          info =  6;
    اگر (n < 0)              info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (order == CblasRowMajor) {
    اگر (Uplo == CblasUpper)         uplo  = 3;
    اگر (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    اگر (incy == 0)          info =  9;
    اگر (incx == 0)          info =  6;
    اگر (n < 0)              info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#endif

  اگر (n == 0) بازگشت;

  اگر ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), پوچ, 0, پوچ, 0);

  اگر ((alpha_r == ZERO) && (alpha_i == ZERO)) بازگشت;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  اگر (incx < 0 ) x -= (n - 1) * incx * 2;
  اگر (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  اگر (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#ifdef SMP
  nthreads = رقم_cpu_avail(2);

  اگر (n < 50) nthreads = 1;
  اگر (nthreads > 2 && n < 500) nthreads = 2;

  اگر (nthreads == 1) {
#endif

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } وگرنه {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  بازگشت;
}
