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
#define ERROR_NAME "QSYR2 "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYR2 "
#បើមិនដូច្នេះ
#define ERROR_NAME "SSYR2 "
#endif

static ចំនួនគត់ (*syr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qsyr2_U, qsyr2_L,
#elif defined(DOUBLE)
  dsyr2_U, dsyr2_L,
#បើមិនដូច្នេះ
  ssyr2_U, ssyr2_L,
#endif
};

#ifdef SMP
static ចំនួនគត់ (*syr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ចំនួនគត់) = {
#ifdef XDOUBLE
  qsyr2_thread_U, qsyr2_thread_L,
#elif defined(DOUBLE)
  dsyr2_thread_U, dsyr2_thread_L,
#បើមិនដូច្នេះ
  ssyr2_thread_U, ssyr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  ចំនួនគត់ uplo;
  FLOAT *buffer;
#ifdef SMP
  ចំនួនគត់ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ប្រសិនបើ (uplo_arg  == 'U') uplo  = 0;
  ប្រសិនបើ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ប្រសិនបើ (lda  < MAX(1, n))   info =  9;
  ប្រសិនបើ (incy == 0)          info =  7;
  ប្រសិនបើ (incx == 0)          info =  5;
  ប្រសិនបើ (n < 0)              info =  2;
  ប្រសិនបើ (uplo  < 0)          info =  1;

  ប្រសិនបើ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

#បើមិនដូច្នេះ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy, FLOAT *a, blasint lda) {

  FLOAT *buffer;
  ចំនួនគត់ uplo;
  blasint info;
#ifdef SMP
  ចំនួនគត់ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ប្រសិនបើ (order == CblasColMajor) {

    ប្រសិនបើ (Uplo == CblasUpper) uplo  = 0;
    ប្រសិនបើ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ប្រសិនបើ (lda  < MAX(1, n))   info =  9;
    ប្រសិនបើ (incy == 0)          info =  7;
    ប្រសិនបើ (incx == 0)          info =  5;
    ប្រសិនបើ (n < 0)              info =  2;
    ប្រសិនបើ (uplo  < 0)          info =  1;
  }

  ប្រសិនបើ (order == CblasRowMajor) {

    ប្រសិនបើ (Uplo == CblasUpper) uplo  = 1;
    ប្រសិនបើ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    ប្រសិនបើ (lda  < MAX(1, n))   info =  9;
    ប្រសិនបើ (incy == 0)          info =  7;
    ប្រសិនបើ (incx == 0)          info =  5;
    ប្រសិនបើ (n < 0)              info =  2;
    ប្រសិនបើ (uplo  < 0)          info =  1;
  }

  ប្រសិនបើ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

#endif

  ប្រសិនបើ (n == 0) ត្រឡប់;

  ប្រសិនបើ (alpha == ZERO) ត្រឡប់;

  IDEBUG_START;

  ប្រសិនបើ (incx == 1 && incy == 1 && n < 100) {
    blasint i;
    ប្រសិនបើ (!uplo) {
      សម្រាប់ (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, NULL, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, NULL, 0);
        a += lda;
      }
    } បើមិនដូច្នេះ {
      សម្រាប់ (i = 0; i < n; i++){
        AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, NULL, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, NULL, 0);
        a += 1 + lda;
      }
    }
    ត្រឡប់;
  }

	  
  FUNCTION_PROFILE_START();

  ប្រសិនបើ (incx < 0 ) x -= (n - 1) * incx;
  ប្រសិនបើ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ប្រសិនបើ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ប្រសិនបើ (n < 50) nthreads = 1;
  ប្រសិនបើ (nthreads > 2 && n < 500) nthreads = 2;

  ប្រសិនបើ (nthreads == 1) {
#endif

    (syr2[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer);

#ifdef SMP
  } បើមិនដូច្នេះ {

    (syr2_thread[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  ត្រឡប់;
}
