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
#không_thì
#define ERROR_NAME "SSYR  "
#endif

static số_nguyên (*syr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qsyr_U, qsyr_L,
#elif defined(DOUBLE)
  dsyr_U, dsyr_L,
#không_thì
  ssyr_U, ssyr_L,
#endif
};

#ifdef SMP
static số_nguyên (*syr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, số_nguyên) = {
#ifdef XDOUBLE
  qsyr_thread_U, qsyr_thread_L,
#elif defined(DOUBLE)
  dsyr_thread_U, dsyr_thread_L,
#không_thì
  ssyr_thread_U, ssyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  char uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  số_nguyên uplo;
  FLOAT *buffer;
#ifdef SMP
  số_nguyên nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  nếu (uplo_arg  == 'U') uplo  = 0;
  nếu (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  nếu (lda  < MAX(1, n))   info =  7;
  nếu (incx == 0)          info =  5;
  nếu (n < 0)              info =  2;
  nếu (uplo  < 0)          info =  1;

  nếu (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

void CNAME(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *a, blasint lda) {

  FLOAT *buffer;
  số_nguyên uplo;
  blasint info;
#ifdef SMP
  số_nguyên nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  nếu (order == CblasColMajor) {

    nếu (Uplo == CblasUpper) uplo  = 0;
    nếu (Uplo == CblasLower) uplo  = 1;

    info = -1;

    nếu (lda  < MAX(1, n))   info =  7;
    nếu (incx == 0)          info =  5;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;

  }

  nếu (order == CblasRowMajor) {

    nếu (Uplo == CblasUpper) uplo  = 1;
    nếu (Uplo == CblasLower) uplo  = 0;

    info = -1;

    nếu (lda  < MAX(1, n))   info =  7;
    nếu (incx == 0)          info =  5;
    nếu (n < 0)              info =  2;
    nếu (uplo  < 0)          info =  1;
  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif

  nếu (n == 0) trả_về;

  nếu (alpha == ZERO) trả_về;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
#nếu 1
  nếu (incx == 1 && n < 100) {
    BLASLONG i;

    nếu (uplo == 0) {
      cho (i = 0; i < n; i++){
        nếu (x[i] != ZERO) {
          AXPYU_K(i + 1, 0, 0, alpha * x[i], x,     1, a, 1, NULL, 0);
        }
        a += lda;
      }  
    } không_thì {
      cho (i = 0; i < n; i++){
        nếu (x[i] != ZERO) {
          AXPYU_K(n - i, 0, 0, alpha * x[i], x + i, 1, a, 1, NULL, 0);
        }
        a += 1 + lda;
      }
    }
    trả_về;
  } 
#endif
  nếu (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  nếu (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  nếu (n < 50) nthreads = 1;
  nếu (nthreads > 2 && n < 500) nthreads = 2;

  nếu (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha, x, incx, a, lda, buffer);

#ifdef SMP
  } không_thì {

    (syr_thread[uplo])(n, alpha, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  trả_về;
}
