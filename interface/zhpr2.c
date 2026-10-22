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
#define ERROR_NAME "XHPR2 "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHPR2 "
#không_thì
#define ERROR_NAME "CHPR2 "
#endif

static số_nguyên (*hpr2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr2_U, xhpr2_L, xhpr2_V, xhpr2_M,
#elif defined(DOUBLE)
  zhpr2_U, zhpr2_L, zhpr2_V, zhpr2_M,
#không_thì
  chpr2_U, chpr2_L, chpr2_V, chpr2_M,
#endif
};

#ifdef SMP
static số_nguyên (*hpr2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, số_nguyên) = {
#ifdef XDOUBLE
  xhpr2_thread_U, xhpr2_thread_L, xhpr2_thread_V, xhpr2_thread_M,
#elif defined(DOUBLE)
  zhpr2_thread_U, zhpr2_thread_L, zhpr2_thread_V, zhpr2_thread_M,
#không_thì
  chpr2_thread_U, chpr2_thread_L, chpr2_thread_V, chpr2_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  char uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
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

  nếu (incy == 0)          info =  7;
  nếu (incx == 0)          info =  5;
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
	   void  *vx, blasint incx,
	   void  *vy, blasint incy,
	   void  *va) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;
  FLOAT* a = (FLOAT*) va;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
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

    nếu (incy == 0)          info =  7;
    nếu (incx == 0)          info =  5;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;
  }

  nếu (order == CblasRowMajor) {
    nếu (Uplo == CblasUpper)         uplo  = 3;
    nếu (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    nếu (incx == 0)          info =  7;
    nếu (incy == 0)          info =  5;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;
  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif

  nếu (n == 0) trả_về;

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
  nếu (n < 100)
	nthreads = 1;
  không_thì
    nthreads = num_cpu_avail(2);

  nếu (nthreads == 1) {
#endif

    (hpr2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, buffer);

#ifdef SMP
  } không_thì {

    (hpr2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  trả_về;
}
