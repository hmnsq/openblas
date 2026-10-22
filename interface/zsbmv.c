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

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "XSBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSBMV "
#否则
#define ERROR_NAME "CSBMV "
#endif

static  整数 (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  xsbmv_U, xsbmv_L,
#elif defined(DOUBLE)
  zsbmv_U, zsbmv_L,
#否则
  csbmv_U, csbmv_L,
#endif
};

#ifdef SMPTEST
static  整数 (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, 整数) = {
#ifdef XDOUBLE
  xsbmv_thread_U, xsbmv_thread_L,
#elif defined(DOUBLE)
  zsbmv_thread_U, zsbmv_thread_L,
#否则
  csbmv_thread_U, csbmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
void NAME(char *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  char uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  整数 uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  整数 nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  如果 (uplo_arg  == 'U') uplo  = 0;
  如果 (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  如果 (incy == 0)          info = 11;
  如果 (incx == 0)          info =  8;
  如果 (lda  < k + 1)       info =  6;
  如果 (k < 0)              info =  3;
  如果 (n < 0)              info =  2;
  如果 (uplo  < 0)          info =  1;

  如果 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    返回;
  }

  如果 (n == 0) 返回;

  如果 ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, blasabs(incy), NULL, 0, NULL, 0);

  如果 ((alpha_r == ZERO) && (alpha_i == ZERO)) 返回;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  如果 (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  如果 (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  如果 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    返回;
  }

#ifdef SMPTEST
  nthreads = num_cpu_avail(2);

  如果 (nthreads == 1) {
#endif

  (sbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, b, incx, c, incy, buffer);

#ifdef SMPTEST
  } 否则 {

    (sbmv_thread[uplo])(n, k, ALPHA, a, lda, b, incx, c, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  返回;
}
