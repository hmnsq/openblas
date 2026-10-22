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
#define ERROR_NAME "XHPR  "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHPR  "
#નહીંતર
#define ERROR_NAME "CHPR  "
#endif

સ્થિર પૂર્ણાંક (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#અથવાજો defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#નહીંતર
  chpr_U, chpr_L, chpr_V, chpr_M,
#endif
};

#ifdef SMP
સ્થિર પૂર્ણાંક (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#અથવાજો defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#નહીંતર
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  અક્ષર uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

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

  જો (incx == 0)          info =  5;
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
	   FLOAT alpha,
	   શૂન્ય  *vx, blasint incx,
	   શૂન્ય  *va) {

  FLOAT* x = (FLOAT*) vx;
  FLOAT* a = (FLOAT*) va;

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

    જો (incx == 0)          info =  5;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Uplo == CblasUpper)         uplo  = 3;
    જો (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    જો (incx == 0)          info =  5;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  જો (n == 0) પાછા;

  જો (alpha == ZERO) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#ifdef SMP
  જો (n < 100)
    nthreads = 1;
  નહીંતર
    nthreads = અંક_cpu_avail(2);

  જો (nthreads == 1) {
#endif

  (hpr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP

  } નહીંતર {

    (hpr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  પાછા;
}
