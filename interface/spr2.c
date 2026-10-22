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

#ifdef XDOUBLE
#define ERROR_NAME "QSPR2 "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DSPR2 "
#இல்லையெனில்
#define ERROR_NAME "SSPR2 "
#endif

நிலையான முழு (*spr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr2_U, qspr2_L,
#அல்லதெனில் defined(DOUBLE)
  dspr2_U, dspr2_L,
#இல்லையெனில்
  sspr2_U, sspr2_L,
#endif
};

#ifdef SMP
நிலையான முழு (*spr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, முழு) = {
#ifdef XDOUBLE
  qspr2_thread_U, qspr2_thread_L,
#அல்லதெனில் defined(DOUBLE)
  dspr2_thread_U, dspr2_thread_L,
#இல்லையெனில்
  sspr2_thread_U, sspr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  எழுத்து uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  முழு uplo;
  FLOAT *buffer;
#ifdef SMP
  முழு nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  எனில் (incy == 0)          info =  7;
  எனில் (incx == 0)          info =  5;
  எனில் (n < 0)              info =  2;
  எனில் (uplo  < 0)          info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order,
	   எண்ணிக்கை CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a) {

  FLOAT *buffer;
  முழு uplo;
  blasint info;
#ifdef SMP
  முழு nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 0;
    எனில் (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    எனில் (incy == 0)          info =  7;
    எனில் (incx == 0)          info =  5;
    எனில் (n < 0)              info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 1;
    எனில் (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    எனில் (incy == 0)          info =  7;
    எனில் (incx == 0)          info =  5;
    எனில் (n < 0)              info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif

  எனில் (n == 0) திரும்பு;

  எனில் (alpha == ZERO) திரும்பு;

  எனில் (incx == 1 && incy == 1 && n < 50) {
    blasint i;
    எனில் (!uplo) {
      க்கு (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, வெற்று, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, வெற்று, 0);
        a += i + 1;
      }
    } இல்லையெனில் {
      க்கு (i = 0; i < n; i++){
	AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, வெற்று, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, வெற்று, 0);
        a += n - i;
      }
    }
    திரும்பு;
  }

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

#ifdef SMP
  nthreads = இலக்கம்_cpu_avail(2);

  எனில் (n < 50) nthreads = 1;
  எனில் (nthreads > 2 && n < 500) nthreads = 2;

  எனில் (nthreads == 1) {
#endif

  (spr2[uplo])(n, alpha, x, incx, y, incy, a, buffer);

#ifdef SMP
  } இல்லையெனில் {

    (spr2_thread[uplo])(n, alpha, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  திரும்பு;
}
