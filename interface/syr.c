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
#define ERROR_NAME "QSYR  "
#elif defined(DOUBLE)
#define ERROR_NAME "DSYR  "
#ᬮᬶᬬ
#define ERROR_NAME "SSYR  "
#endif

static puur_naa_nka (*syr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qsyr_U, qsyr_L,
#elif defined(DOUBLE)
  dsyr_U, dsyr_L,
#ᬮᬶᬬ
  ssyr_U, ssyr_L,
#endif
};

#ifdef SMP
static puur_naa_nka (*syr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  qsyr_thread_U, qsyr_thread_L,
#elif defined(DOUBLE)
  dsyr_thread_U, dsyr_thread_L,
#ᬮᬶᬬ
  ssyr_thread_U, ssyr_thread_L,
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
  puur_naa_nka uplo;
  FLOAT *buffer;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ᬬᬾᬦ᭄ (uplo_arg  == 'U') uplo  = 0;
  ᬬᬾᬦ᭄ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ᬬᬾᬦ᭄ (lda  < MAX(1, n))   info =  7;
  ᬬᬾᬦ᭄ (incx == 0)          info =  5;
  ᬬᬾᬦ᭄ (n < 0)              info =  2;
  ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;

  ᬬᬾᬦ᭄ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#ᬮᬶᬬ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *a, blasint lda) {

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ᬬᬾᬦ᭄ (order == CblasColMajor) {

    ᬬᬾᬦ᭄ (Uplo == CblasUpper) uplo  = 0;
    ᬬᬾᬦ᭄ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ᬬᬾᬦ᭄ (lda  < MAX(1, n))   info =  7;
    ᬬᬾᬦ᭄ (incx == 0)          info =  5;
    ᬬᬾᬦ᭄ (n < 0)              info =  2;
    ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;

  }

  ᬬᬾᬦ᭄ (order == CblasRowMajor) {

    ᬬᬾᬦ᭄ (Uplo == CblasUpper) uplo  = 1;
    ᬬᬾᬦ᭄ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    ᬬᬾᬦ᭄ (lda  < MAX(1, n))   info =  7;
    ᬬᬾᬦ᭄ (incx == 0)          info =  5;
    ᬬᬾᬦ᭄ (n < 0)              info =  2;
    ᬬᬾᬦ᭄ (uplo  < 0)          info =  1;
  }

  ᬬᬾᬦ᭄ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#endif

  ᬬᬾᬦ᭄ (n == 0) ᬯᬮᬶ;

  ᬬᬾᬦ᭄ (alpha == ZERO) ᬯᬮᬶ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
#ᬬᬾᬦ᭄ 1
  ᬬᬾᬦ᭄ (incx == 1 && n < 100) {
    BLASLONG i;

    ᬬᬾᬦ᭄ (uplo == 0) {
      for (i = 0; i < n; i++){
        ᬬᬾᬦ᭄ (x[i] != ZERO) {
          AXPYU_K(i + 1, 0, 0, alpha * x[i], x,     1, a, 1, NULL, 0);
        }
        a += lda;
      }  
    } ᬮᬶᬬ {
      for (i = 0; i < n; i++){
        ᬬᬾᬦ᭄ (x[i] != ZERO) {
          AXPYU_K(n - i, 0, 0, alpha * x[i], x + i, 1, a, 1, NULL, 0);
        }
        a += 1 + lda;
      }
    }
    ᬯᬮᬶ;
  } 
#endif
  ᬬᬾᬦ᭄ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᬬᬾᬦ᭄ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᬯᬮᬶ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ᬬᬾᬦ᭄ (n < 50) nthreads = 1;
  ᬬᬾᬦ᭄ (nthreads > 2 && n < 500) nthreads = 2;

  ᬬᬾᬦ᭄ (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha, x, incx, a, lda, buffer);

#ifdef SMP
  } ᬮᬶᬬ {

    (syr_thread[uplo])(n, alpha, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ᬯᬮᬶ;
}
