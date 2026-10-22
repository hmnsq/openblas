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
#define ERROR_NAME "QSPR  "
#elif defined(DOUBLE)
#define ERROR_NAME "DSPR  "
#ᯀᯮᯝᯂᯪ
#define ERROR_NAME "SSPR  "
#endif

static puur_naa_nka (*spr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr_U, qspr_L,
#elif defined(DOUBLE)
  dspr_U, dspr_L,
#ᯀᯮᯝᯂᯪ
  sspr_U, sspr_L,
#endif
};

#ifdef SMP
static puur_naa_nka (*spr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, puur_naa_nka) = {
#ifdef XDOUBLE
  qspr_thread_U, qspr_thread_L,
#elif defined(DOUBLE)
  dspr_thread_U, dspr_thread_L,
#ᯀᯮᯝᯂᯪ
  sspr_thread_U, sspr_thread_L,
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

  ᯀᯮᯝ (incx == 0)          info =  5;
  ᯀᯮᯝ (n < 0)              info =  2;
  ᯀᯮᯝ (uplo  < 0)          info =  1;

  ᯀᯮᯝ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᯀᯮᯞ;
  }

#ᯀᯮᯝᯂᯪ

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *a) {

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  ᯀᯮᯝ (order == CblasColMajor) {
    ᯀᯮᯝ (Uplo == CblasUpper)         uplo  = 0;
    ᯀᯮᯝ (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    ᯀᯮᯝ (incx == 0)          info =  5;
    ᯀᯮᯝ (n < 0)              info =  2;
    ᯀᯮᯝ (uplo  < 0)          info =  1;
  }

  ᯀᯮᯝ (order == CblasRowMajor) {
    ᯀᯮᯝ (Uplo == CblasUpper)         uplo  = 1;
    ᯀᯮᯝ (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    ᯀᯮᯝ (incx == 0)          info =  5;
    ᯀᯮᯝ (n < 0)              info =  2;
    ᯀᯮᯝ (uplo  < 0)          info =  1;
  }

  ᯀᯮᯝ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᯀᯮᯞ;
  }

#endif

  ᯀᯮᯝ (n == 0) ᯀᯮᯞ;

  ᯀᯮᯝ (alpha == ZERO) ᯀᯮᯞ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᯀᯮᯝ (incx == 1 && n <100) {
    blasint i;
    ᯀᯮᯝ (uplo==0) {
      for (i = 0; i < n; i++){
        ᯀᯮᯝ (x[i] != ZERO) {
          AXPYU_K(i + 1, 0, 0, alpha * x[i], x,     1, a, 1, NULL, 0);
        }
        a += i + 1;
      }
    } ᯀᯮᯝᯂᯪ { 
      for (i = 0; i < n; i++){
        ᯀᯮᯝ (x[i] != ZERO) {
          AXPYU_K(n - i, 0, 0, alpha * x[i], x + i, 1, a, 1, NULL, 0);
        }
        a += n - i;
      }
    }
    ᯀᯮᯞ;
  }

  ᯀᯮᯝ (incx < 0 ) x -= (n - 1) * incx;

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

  (spr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP
  } ᯀᯮᯝᯂᯪ {

    (spr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ᯀᯮᯞ;
}
