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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <ctype.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHPMV "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHPMV "
#નહીંતર
#define ERROR_NAME "CHPMV "
#endif

સ્થિર  પૂર્ણાંક (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, શૂન્ય *) = {
#ifdef XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#અથવાજો defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#નહીંતર
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#endif
};

#ifdef SMP
સ્થિર  પૂર્ણાંક (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#અથવાજો defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#નહીંતર
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  અક્ષર uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  પૂર્ણાંક uplo;
  FLOAT *buffer;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  જો (incy == 0)          info =  9;
  જો (incx == 0)          info =  6;
  જો (n < 0)              info =  2;
  જો (uplo  < 0)          info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   ગણના CBLAS_UPLO Uplo,
	   blasint n,
	   શૂન્ય *VALPHA,
	   શૂન્ય  *va,
	   શૂન્ય  *vx, blasint incx,
	   શૂન્ય *VBETA,
	   શૂન્ય  *vy, blasint incy){

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
  પૂર્ણાંક uplo;
  blasint info;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  જો (order == CblasColMajor) {
    જો (Uplo == CblasUpper)         uplo  = 0;
    જો (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    જો (incy == 0)          info =  9;
    જો (incx == 0)          info =  6;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Uplo == CblasUpper)         uplo  = 3;
    જો (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    જો (incy == 0)          info =  9;
    જો (incx == 0)          info =  6;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  જો (n == 0) પાછા;

  જો ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), ખાલી, 0, ખાલી, 0);

  જો ((alpha_r == ZERO) && (alpha_i == ZERO)) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0 ) x -= (n - 1) * incx * 2;
  જો (incy < 0 ) y -= (n - 1) * incy * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#ifdef SMP
  nthreads = અંક_cpu_avail(2);

  જો (n < 50) nthreads = 1;
  જો (nthreads > 2 && n < 500) nthreads = 2;

  જો (nthreads == 1) {
#endif

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, buffer);

#ifdef SMP
  } નહીંતર {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  પાછા;
}
