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
#define ERROR_NAME "QSYMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYMV "
#དེ་མིན
#define ERROR_NAME "SSYMV "
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na uplo_arg = *UPLO;
  blasint n	= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  ཧྲིལ (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ifdef SMP
  ཧྲིལ (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ཧྲིལ) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#endif

  blasint info;
  ཧྲིལ uplo;
  FLOAT *buffer;
#ifdef SMP
  ཧྲིལ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  གལ་ཏེ (uplo_arg  == 'U') uplo  = 0;
  གལ་ཏེ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  གལ་ཏེ (incy == 0)          info = 10;
  གལ་ཏེ (incx == 0)          info =  7;
  གལ་ཏེ (lda  < MAX(1, n))   info =  5;
  གལ་ཏེ (n < 0)              info =  2;
  གལ་ཏེ (uplo  < 0)          info =  1;

  གལ་ཏེ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

#དེ་མིན

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha,
	   FLOAT *a, blasint lda, FLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy) {

  FLOAT *buffer;
  ཧྲིལ uplo;
  blasint info;
#ifdef SMP
  ཧྲིལ nthreads;
#endif

  ཧྲིལ (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ifdef SMP
  ཧྲིལ (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ཧྲིལ) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  གལ་ཏེ (order == CblasColMajor) {

    གལ་ཏེ (Uplo == CblasUpper) uplo  = 0;
    གལ་ཏེ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    གལ་ཏེ (incy == 0)          info = 10;
    གལ་ཏེ (incx == 0)          info =  7;
    གལ་ཏེ (lda  < MAX(1, n))   info =  5;
    གལ་ཏེ (n < 0)              info =  2;
    གལ་ཏེ (uplo  < 0)          info =  1;
  }

  གལ་ཏེ (order == CblasRowMajor) {

    གལ་ཏེ (Uplo == CblasUpper) uplo  = 1;
    གལ་ཏེ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    གལ་ཏེ (incy == 0)          info = 10;
    གལ་ཏེ (incx == 0)          info =  7;
    གལ་ཏེ (lda  < MAX(1, n))   info =  5;
    གལ་ཏེ (n < 0)              info =  2;
    གལ་ཏེ (uplo  < 0)          info =  1;
  }

  གལ་ཏེ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

#endif

  གལ་ཏེ (n == 0) ཕྱིར;

  གལ་ཏེ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  གལ་ཏེ (alpha == ZERO) ཕྱིར;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  གལ་ཏེ (incx < 0 ) x -= (n - 1) * incx;
  གལ་ཏེ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  གལ་ཏེ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ཕྱིར;
  }

#ifdef SMP
  གལ་ཏེ (n <200)
	  nthreads=1;
  དེ་མིན
  nthreads = num_cpu_avail(2);

  གལ་ཏེ (nthreads == 1) {
#endif

  (symv[uplo])(n, n, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } དེ་མིན {

    (symv_thread[uplo])(n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n,  2 * n * n);

  IDEBUG_END;

  ཕྱིར;
}
