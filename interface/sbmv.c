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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <ctype.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "QSBMV "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DSBMV "
#இல்லையெனில்
#define ERROR_NAME "SSBMV "
#endif

நிலையான  முழு (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, வெறுமை *) = {
#ifdef XDOUBLE
  qsbmv_U, qsbmv_L,
#அல்லதெனில் defined(DOUBLE)
  dsbmv_U, dsbmv_L,
#இல்லையெனில்
  ssbmv_U, ssbmv_L,
#endif
};

#ifdef SMPTEST
நிலையான  முழு (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, முழு) = {
#ifdef XDOUBLE
  qsbmv_thread_U, qsbmv_thread_L,
#அல்லதெனில் defined(DOUBLE)
  dsbmv_thread_U, dsbmv_thread_L,
#இல்லையெனில்
  ssbmv_thread_U, ssbmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  எழுத்து uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  முழு uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  முழு nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  எனில் (incy == 0)          info = 11;
  எனில் (incx == 0)          info =  8;
  எனில் (lda  < k + 1)       info =  6;
  எனில் (k < 0)              info =  3;
  எனில் (n < 0)              info =  2;
  எனில் (uplo  < 0)          info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   எண்ணிக்கை CBLAS_UPLO Uplo,
	   blasint n, blasint k,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  முழு uplo;
  blasint info;
#ifdef SMPTEST
  முழு nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 0;
    எனில் (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    எனில் (incy == 0)          info = 11;
    எனில் (incx == 0)          info =  8;
    எனில் (lda  < k + 1)       info =  6;
    எனில் (k < 0)              info =  3;
    எனில் (n < 0)              info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 1;
    எனில் (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    எனில் (incy == 0)          info = 11;
    எனில் (incx == 0)          info =  8;
    எனில் (lda  < k + 1)       info =  6;
    எனில் (k < 0)              info =  3;
    எனில் (n < 0)              info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif

  எனில் (n == 0) திரும்பு;

  எனில் (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), வெற்று, 0, வெற்று, 0);

  எனில் (alpha == ZERO) திரும்பு;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  எனில் (incx < 0 ) x -= (n - 1) * incx;
  எனில் (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#ifdef SMPTEST
  nthreads = இலக்கம்_cpu_avail(2);

  எனில் (nthreads == 1) {
#endif

  (sbmv[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMPTEST
  } இல்லையெனில் {

    (sbmv_thread[uplo])(n, k, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * k / 2 + n, n * k);

  IDEBUG_END;

  திரும்பு;
}
