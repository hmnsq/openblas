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
#define ERROR_NAME "ZSYR2 "
#ᯀᯮᯝᯂᯪ
#define ERROR_NAME "CSYR2 "
#endif

static puur_naa_nka (*syr2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr2_U, xsyr2_L,
#elif defined(DOUBLE)
  zsyr2_U, zsyr2_L,
#ᯀᯮᯝᯂᯪ
  csyr2_U, csyr2_L,
#endif
};

#ifdef SMP
static puur_naa_nka (*syr2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  xsyr2_thread_U, xsyr2_thread_L,
#elif defined(DOUBLE)
  zsyr2_thread_U, zsyr2_thread_L,
#ᯀᯮᯝᯂᯪ
  csyr2_thread_U, csyr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
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

  ᯀᯮᯝ (uplo_arg  == 'U') uplo  = 0;
  ᯀᯮᯝ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᯀᯮᯝ (lda  < MAX(1, n))   info =  9;
  ᯀᯮᯝ (incy == 0)          info =  7;
  ᯀᯮᯝ (incx == 0)          info =  5;
  ᯀᯮᯝ (n < 0)              info =  2;
  ᯀᯮᯝ (uplo  < 0)          info =  1;

  ᯀᯮᯝ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᯀᯮᯞ;
  }

  ᯀᯮᯝ (n == 0) ᯀᯮᯞ;

  ᯀᯮᯝ ((alpha_r == ZERO) && (alpha_i == ZERO)) ᯀᯮᯞ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᯀᯮᯝ (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;
  ᯀᯮᯝ (incy < 0 ) y -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᯀᯮᯝ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᯀᯮᯞ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ᯀᯮᯝ (n < 50) nthreads = 1;
  ᯀᯮᯝ (nthreads > 2 && n < 500) nthreads = 2;

  ᯀᯮᯝ (nthreads == 1) {
#endif

    (syr2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, lda, buffer);

#ifdef SMP
  } ᯀᯮᯝᯂᯪ {

    (syr2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  ᯀᯮᯞ;
}
