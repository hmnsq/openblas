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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHPMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#không_thì
#define ERROR_NAME "CHPMV "
#endif

static  số_nguyên (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#elif defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#không_thì
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
static  số_nguyên (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, số_nguyên) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#elif defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#không_thì
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  char uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  số_nguyên uplo;
  FLOAT *buffer;
#ifdef SMP
  số_nguyên nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  nếu (uplo_arg  == 'U') uplo  = 0;
  nếu (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  nếu (incy == 0)          info =  9;
  nếu (incx == 0)          info =  6;
  nếu (n < 0)              info =  2;
  nếu (uplo  < 0)          info =  1;

  nếu (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

void CNAME(enum CBLAS_ORDER order,
	   enum CBLAS_UPLO Uplo,
	   blasint n,
	   void *VALPHA,
	   void  *va,
	   void  *vx, blasint incx,
	   void *VBETA,
	   void  *vy, blasint incy){

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
  số_nguyên uplo;
  blasint info;
#ifdef SMP
  số_nguyên nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  nếu (order == CblasColMajor) {
    nếu (Uplo == CblasUpper)         uplo  = 0;
    nếu (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    nếu (incy == 0)          info =  9;
    nếu (incx == 0)          info =  6;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;
  }

  nếu (order == CblasRowMajor) {
    nếu (Uplo == CblasUpper)         uplo  = 3;
    nếu (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    nếu (incy == 0)          info =  9;
    nếu (incx == 0)          info =  6;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;
  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif

  nếu (n == 0) trả_về;

  nếu ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), NULL, 0, NULL, 0);

  nếu ((alpha_r == ZERO) && (alpha_i == ZERO)) trả_về;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  nếu (incx < 0 ) x -= (n - 1) * incx * 2;
  nếu (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  nếu (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  nếu (n < 50) nthreads = 1;
  nếu (nthreads > 2 && n < 500) nthreads = 2;

  nếu (nthreads == 1) {
#endif

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } không_thì {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  trả_về;
}
