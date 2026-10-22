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
#define ERROR_NAME "XHPR  "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHPR  "
#มิฉะนั้น
#define ERROR_NAME "CHPR  "
#endif

static จำนวนเต็ม (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#elif defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#มิฉะนั้น
  chpr_U, chpr_L, chpr_V, chpr_M,
#endif
};

#ifdef SMP
static จำนวนเต็ม (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, จำนวนเต็ม) = {
#ifdef XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#elif defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#มิฉะนั้น
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

  blasint info;
  จำนวนเต็ม uplo;
  FLOAT *buffer;
#ifdef SMP
  จำนวนเต็ม nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ถ้า (uplo_arg  == 'U') uplo  = 0;
  ถ้า (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ถ้า (incx == 0)          info =  5;
  ถ้า (n < 0)              info =  2;
  ถ้า (uplo  < 0)          info =  1;

  ถ้า (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#มิฉะนั้น

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   void  *vx, blasint incx,
	   void  *va) {

  FLOAT* x = (FLOAT*) vx;
  FLOAT* a = (FLOAT*) va;

  FLOAT *buffer;
  จำนวนเต็ม uplo;
  blasint info;
#ifdef SMP
  จำนวนเต็ม nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ถ้า (order == CblasColMajor) {
    ถ้า (Uplo == CblasUpper)         uplo  = 0;
    ถ้า (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ถ้า (incx == 0)          info =  5;
    ถ้า (n < 0)              info =  2;
    ถ้า (uplo  < 0)          info =  1;
  }

  ถ้า (order == CblasRowMajor) {
    ถ้า (Uplo == CblasUpper)         uplo  = 3;
    ถ้า (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    ถ้า (incx == 0)          info =  5;
    ถ้า (n < 0)              info =  2;
    ถ้า (uplo  < 0)          info =  1;
  }

  ถ้า (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#endif

  ถ้า (n == 0) คืน;

  ถ้า (alpha == ZERO) คืน;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ถ้า (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ถ้า (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#ifdef SMP
  ถ้า (n < 100)
    nthreads = 1;
  มิฉะนั้น
    nthreads = num_cpu_avail(2);

  ถ้า (nthreads == 1) {
#endif

  (hpr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP

  } มิฉะนั้น {

    (hpr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  คืน;
}
