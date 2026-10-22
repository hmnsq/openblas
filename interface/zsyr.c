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
#define ERROR_NAME "XSYR  "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR  "
#ᤀᤡ᤺ᤆᤡ
#define ERROR_NAME "CSYR  "
#endif

static puur_naa_nka (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr_U, xsyr_L,
#elif defined(DOUBLE)
  zsyr_U, zsyr_L,
#ᤀᤡ᤺ᤆᤡ
  csyr_U, csyr_L,
#endif
};

#ifdef SMP
static puur_naa_nka (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#elif defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#ᤀᤡ᤺ᤆᤡ
  csyr_thread_U, csyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  puur_naa_nka uplo;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ᤀᤡ᤺ (uplo_arg  == 'U') uplo  = 0;
  ᤀᤡ᤺ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᤀᤡ᤺ (lda  < MAX(1, n))   info =  7;
  ᤀᤡ᤺ (incx == 0)          info =  5;
  ᤀᤡ᤺ (n < 0)              info =  2;
  ᤀᤡ᤺ (uplo  < 0)          info =  1;

  ᤀᤡ᤺ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }


#ᤀᤡ᤺ᤆᤡ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, puur_naa_nka n, void* valpha, FLOAT *x, puur_naa_nka incx, FLOAT *a, puur_naa_nka lda) {

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ᤀᤡ᤺ (order == CblasColMajor) {

    ᤀᤡ᤺ (Uplo == CblasUpper) uplo  = 0;
    ᤀᤡ᤺ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  7;
    ᤀᤡ᤺ (incx == 0)          info =  5;
    ᤀᤡ᤺ (n < 0)              info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;

  }

  ᤀᤡ᤺ (order == CblasRowMajor) {

    ᤀᤡ᤺ (Uplo == CblasUpper) uplo  = 1;
    ᤀᤡ᤺ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  7;
    ᤀᤡ᤺ (incx == 0)          info =  5;
    ᤀᤡ᤺ (n < 0)              info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;
  }

  ᤀᤡ᤺ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#endif

  ᤀᤡ᤺ (n == 0) ᤀᤡᤴ;

  ᤀᤡ᤺ ((alpha_r == ZERO) && (alpha_i == ZERO)) ᤀᤡᤴ;

  ᤀᤡ᤺ (incx == 1 && n < 50) {
    blasint i;
    ᤀᤡ᤺ (!uplo) {
      for (i = 0; i < n; i++){
        ᤀᤡ᤺ ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, NULL, 0);
        }
        a += lda * 2;
      }
    } ᤀᤡ᤺ᤆᤡ {
      for (i = 0; i < n; i++){
        ᤀᤡ᤺ ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, NULL, 0);
        }
        a += 2 + lda * 2;
      }
    }
    ᤀᤡᤴ;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᤀᤡ᤺ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᤀᤡ᤺ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ᤀᤡ᤺ (n < 50) nthreads = 1;
  ᤀᤡ᤺ (nthreads > 2 && n < 500) nthreads = 2;

  ᤀᤡ᤺ (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, buffer);

#ifdef SMP
  } ᤀᤡ᤺ᤆᤡ {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ᤀᤡᤴ;
}
