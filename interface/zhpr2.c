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
#ᱵᱟᱝ
#define ERROR_NAME "CHPR2 "
#endif

static puur_naa_nka (*hpr2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr2_U, xhpr2_L, xhpr2_V, xhpr2_M,
#elif defined(DOUBLE)
  zhpr2_U, zhpr2_L, zhpr2_V, zhpr2_M,
#ᱵᱟᱝ
  chpr2_U, chpr2_L, chpr2_V, chpr2_M,
#endif
};

#ifdef SMP
static puur_naa_nka (*hpr2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  xhpr2_thread_U, xhpr2_thread_L, xhpr2_thread_V, xhpr2_thread_M,
#elif defined(DOUBLE)
  zhpr2_thread_U, zhpr2_thread_L, zhpr2_thread_V, zhpr2_thread_M,
#ᱵᱟᱝ
  chpr2_thread_U, chpr2_thread_L, chpr2_thread_V, chpr2_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  puur_naa_nka uplo;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ᱡᱩᱫᱤ (uplo_arg  == 'U') uplo  = 0;
  ᱡᱩᱫᱤ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᱡᱩᱫᱤ (incy == 0)          info =  7;
  ᱡᱩᱫᱤ (incx == 0)          info =  5;
  ᱡᱩᱫᱤ (n < 0)              info =  2;
  ᱡᱩᱫᱤ (uplo  < 0)          info =  1;

  ᱡᱩᱫᱤ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᱨᱩᱣᱟᱹᱲ;
  }

#ᱵᱟᱝ

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_UPLO Uplo,
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
  puur_naa_nka uplo;
  blasint info;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ᱡᱩᱫᱤ (order == CblasColMajor) {
    ᱡᱩᱫᱤ (Uplo == CblasUpper)         uplo  = 0;
    ᱡᱩᱫᱤ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ᱡᱩᱫᱤ (incy == 0)          info =  7;
    ᱡᱩᱫᱤ (incx == 0)          info =  5;
    ᱡᱩᱫᱤ (n < 0)              info =  2;
    ᱡᱩᱫᱤ (uplo  < 0)          info =  1;
  }

  ᱡᱩᱫᱤ (order == CblasRowMajor) {
    ᱡᱩᱫᱤ (Uplo == CblasUpper)         uplo  = 3;
    ᱡᱩᱫᱤ (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ᱡᱩᱫᱤ (incx == 0)          info =  7;
    ᱡᱩᱫᱤ (incy == 0)          info =  5;
    ᱡᱩᱫᱤ (n < 0)              info =  2;
    ᱡᱩᱫᱤ (uplo  < 0)          info =  1;
  }

  ᱡᱩᱫᱤ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᱨᱩᱣᱟᱹᱲ;
  }

#endif

  ᱡᱩᱫᱤ (n == 0) ᱨᱩᱣᱟᱹᱲ;

  ᱡᱩᱫᱤ ((alpha_r == ZERO) && (alpha_i == ZERO)) ᱨᱩᱣᱟᱹᱲ;


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᱡᱩᱫᱤ (incx < 0 ) x -= (n - 1) * incx * 2;
  ᱡᱩᱫᱤ (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᱡᱩᱫᱤ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᱨᱩᱣᱟᱹᱲ;
  }

#ifdef SMP
  ᱡᱩᱫᱤ (n < 100)
	nthreads = 1;
  ᱵᱟᱝ
    nthreads = num_cpu_avail(2);

  ᱡᱩᱫᱤ (nthreads == 1) {
#endif

    (hpr2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, buffer);

#ifdef SMP
  } ᱵᱟᱝ {

    (hpr2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  ᱨᱩᱣᱟᱹᱲ;
}
