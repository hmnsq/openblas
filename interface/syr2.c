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
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "SSYR2 "
#endif

static ຈຳນວນເຕັມ (*syr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  qsyr2_U, qsyr2_L,
#elif defined(DOUBLE)
  dsyr2_U, dsyr2_L,
#ບໍ່ດັ່ງນັ້ນ
  ssyr2_U, ssyr2_L,
#endif
};

#ifdef SMP
static ຈຳນວນເຕັມ (*syr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ຈຳນວນເຕັມ) = {
#ifdef XDOUBLE
  qsyr2_thread_U, qsyr2_thread_L,
#elif defined(DOUBLE)
  dsyr2_thread_U, dsyr2_thread_L,
#ບໍ່ດັ່ງນັ້ນ
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
  ຈຳນວນເຕັມ uplo;
  FLOAT *buffer;
#ifdef SMP
  ຈຳນວນເຕັມ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ຖ້າ (uplo_arg  == 'U') uplo  = 0;
  ຖ້າ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ຖ້າ (lda  < MAX(1, n))   info =  9;
  ຖ້າ (incy == 0)          info =  7;
  ຖ້າ (incx == 0)          info =  5;
  ຖ້າ (n < 0)              info =  2;
  ຖ້າ (uplo  < 0)          info =  1;

  ຖ້າ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#ບໍ່ດັ່ງນັ້ນ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy, FLOAT *a, blasint lda) {

  FLOAT *buffer;
  ຈຳນວນເຕັມ uplo;
  blasint info;
#ifdef SMP
  ຈຳນວນເຕັມ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ຖ້າ (order == CblasColMajor) {

    ຖ້າ (Uplo == CblasUpper) uplo  = 0;
    ຖ້າ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    ຖ້າ (lda  < MAX(1, n))   info =  9;
    ຖ້າ (incy == 0)          info =  7;
    ຖ້າ (incx == 0)          info =  5;
    ຖ້າ (n < 0)              info =  2;
    ຖ້າ (uplo  < 0)          info =  1;
  }

  ຖ້າ (order == CblasRowMajor) {

    ຖ້າ (Uplo == CblasUpper) uplo  = 1;
    ຖ້າ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    ຖ້າ (lda  < MAX(1, n))   info =  9;
    ຖ້າ (incy == 0)          info =  7;
    ຖ້າ (incx == 0)          info =  5;
    ຖ້າ (n < 0)              info =  2;
    ຖ້າ (uplo  < 0)          info =  1;
  }

  ຖ້າ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#endif

  ຖ້າ (n == 0) ສົ່ງຄືນ;

  ຖ້າ (alpha == ZERO) ສົ່ງຄືນ;

  IDEBUG_START;

  ຖ້າ (incx == 1 && incy == 1 && n < 100) {
    blasint i;
    ຖ້າ (!uplo) {
      ສຳລັບ (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, NULL, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, NULL, 0);
        a += lda;
      }
    } ບໍ່ດັ່ງນັ້ນ {
      ສຳລັບ (i = 0; i < n; i++){
        AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, NULL, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, NULL, 0);
        a += 1 + lda;
      }
    }
    ສົ່ງຄືນ;
  }

	  
  FUNCTION_PROFILE_START();

  ຖ້າ (incx < 0 ) x -= (n - 1) * incx;
  ຖ້າ (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ຖ້າ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ຖ້າ (n < 50) nthreads = 1;
  ຖ້າ (nthreads > 2 && n < 500) nthreads = 2;

  ຖ້າ (nthreads == 1) {
#endif

    (syr2[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer);

#ifdef SMP
  } ບໍ່ດັ່ງນັ້ນ {

    (syr2_thread[uplo])(n, alpha, x, incx, y, incy, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  ສົ່ງຄືນ;
}
