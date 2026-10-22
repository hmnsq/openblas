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
#ᤀᤡ᤺ᤆᤡ
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

  puur_naa_nka (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ifdef SMP
  puur_naa_nka (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#endif

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

  ᤀᤡ᤺ (incy == 0)          info = 10;
  ᤀᤡ᤺ (incx == 0)          info =  7;
  ᤀᤡ᤺ (lda  < MAX(1, n))   info =  5;
  ᤀᤡ᤺ (n < 0)              info =  2;
  ᤀᤡ᤺ (uplo  < 0)          info =  1;

  ᤀᤡ᤺ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#ᤀᤡ᤺ᤆᤡ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha,
	   FLOAT *a, blasint lda, FLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy) {

  FLOAT *buffer;
  puur_naa_nka uplo;
  blasint info;
#ifdef SMP
  puur_naa_nka nthreads;
#endif

  puur_naa_nka (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ifdef SMP
  puur_naa_nka (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, puur_naa_nka) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ᤀᤡ᤺ (order == CblasColMajor) {

    ᤀᤡ᤺ (Uplo == CblasUpper) uplo  = 0;
    ᤀᤡ᤺ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ᤀᤡ᤺ (incy == 0)          info = 10;
    ᤀᤡ᤺ (incx == 0)          info =  7;
    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  5;
    ᤀᤡ᤺ (n < 0)              info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;
  }

  ᤀᤡ᤺ (order == CblasRowMajor) {

    ᤀᤡ᤺ (Uplo == CblasUpper) uplo  = 1;
    ᤀᤡ᤺ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    ᤀᤡ᤺ (incy == 0)          info = 10;
    ᤀᤡ᤺ (incx == 0)          info =  7;
    ᤀᤡ᤺ (lda  < MAX(1, n))   info =  5;
    ᤀᤡ᤺ (n < 0)              info =  2;
    ᤀᤡ᤺ (uplo  < 0)          info =  1;
  }

  ᤀᤡ᤺ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#endif

  ᤀᤡ᤺ (n == 0) ᤀᤡᤴ;

  ᤀᤡ᤺ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  ᤀᤡ᤺ (alpha == ZERO) ᤀᤡᤴ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ᤀᤡ᤺ (incx < 0 ) x -= (n - 1) * incx;
  ᤀᤡ᤺ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ᤀᤡ᤺ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ᤀᤡᤴ;
  }

#ifdef SMP
  ᤀᤡ᤺ (n <200)
	  nthreads=1;
  ᤀᤡ᤺ᤆᤡ
  nthreads = num_cpu_avail(2);

  ᤀᤡ᤺ (nthreads == 1) {
#endif

  (symv[uplo])(n, n, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ᤀᤡ᤺ᤆᤡ {

    (symv_thread[uplo])(n, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n,  2 * n * n);

  IDEBUG_END;

  ᤀᤡᤴ;
}
