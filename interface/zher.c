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
#define ERROR_NAME "XHER  "
#elif defined(DOUBLE)
#define ERROR_NAME "ZHER  "
#དེ་མིན
#define ERROR_NAME "CHER  "
#endif

static ཧྲིལ (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#elif defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#དེ་མིན
  cher_U, cher_L, cher_V, cher_M,
#endif
};

#ifdef SMP
static ཧྲིལ (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ཧྲིལ) = {
#ifdef XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#elif defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#དེ་མིན
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  ཧྲིལ uplo;
  FLOAT *buffer;
#ifdef SMP
  ཧྲིལ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  གལ་སྲིད (uplo_arg  == 'U') uplo  = 0;
  གལ་སྲིད (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  གལ་སྲིད (lda  < MAX(1, n))   info =  7;
  གལ་སྲིད (incx == 0)          info =  5;
  གལ་སྲིད (n < 0)              info =  2;
  གལ་སྲིད (uplo  < 0)          info =  1;

  གལ་སྲིད (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

#དེ་མིན

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha, void *vx, blasint incx, void *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *buffer;
  ཧྲིལ uplo;
  blasint info;
#ifdef SMP
  ཧྲིལ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  གལ་སྲིད (order == CblasColMajor) {

    གལ་སྲིད (Uplo == CblasUpper) uplo  = 0;
    གལ་སྲིད (Uplo == CblasLower) uplo  = 1;

    info = -1;

    གལ་སྲིད (lda  < MAX(1, n))   info =  7;
    གལ་སྲིད (incx == 0)          info =  5;
    གལ་སྲིད (n < 0)              info =  2;
    གལ་སྲིད (uplo  < 0)          info =  1;

  }

  གལ་སྲིད (order == CblasRowMajor) {

    གལ་སྲིད (Uplo == CblasUpper) uplo  = 3;
    གལ་སྲིད (Uplo == CblasLower) uplo  = 2;

    info = -1;

    གལ་སྲིད (lda  < MAX(1, n))   info =  7;
    གལ་སྲིད (incx == 0)          info =  5;
    གལ་སྲིད (n < 0)              info =  2;
    གལ་སྲིད (uplo  < 0)          info =  1;
  }

  གལ་སྲིད (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

#endif

  གལ་སྲིད (n == 0) ལོག;

  གལ་སྲིད (alpha == ZERO) ལོག;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  གལ་སྲིད (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  གལ་སྲིད (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ལོག;
  }

#ifdef SMP
  གལ་སྲིད (n < 100)
    nthreads = 1;
  དེ་མིན
    nthreads = num_cpu_avail(2);

  གལ་སྲིད (nthreads == 1) {
#endif

  (her[uplo])(n, alpha, x, incx, a, lda, buffer);

#ifdef SMP
  } དེ་མིན {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ལོག;
}
