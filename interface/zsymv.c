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
#define ERROR_NAME "XSYMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYMV "
#បើមិនដូច្នេះ
#define ERROR_NAME "CSYMV "
#endif

OPENBLAS_EXPORT

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  ចំនួនគត់ (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ifdef SMP
  ចំនួនគត់ (*symv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ចំនួនគត់) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#endif

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

  ប្រសិនបើ (incy == 0)          info = 10;
  ប្រសិនបើ (incx == 0)          info =  7;
  ប្រសិនបើ (lda  < MAX(1, n))   info =  5;
  ប្រសិនបើ (n < 0)              info =  2;
  ប្រសិនបើ (uplo  < 0)          info =  1;

  ប្រសិនបើ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ត្រឡប់;
  }

  ប្រសិនបើ (n == 0) ត្រឡប់;

  ប្រសិនបើ ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, abs(incy), NULL, 0, NULL, 0);

  ប្រសិនបើ ((alpha_r == ZERO) && (alpha_i == ZERO)) ត្រឡប់;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ប្រសិនបើ (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  ប្រសិនបើ (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

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

  (symv[uplo])(n, n, alpha_r, alpha_i, a, lda, b, incx, c, incy, buffer);

#ifdef SMP
  } បើមិនដូច្នេះ {

    (symv_thread[uplo])(n, ALPHA, a, lda, b, incx, c, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n,  2 * n * n);

  IDEBUG_END;

  ត្រឡប់;
}
