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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QSPR2 "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DSPR2 "
#אחרת
#define ERROR_NAME "SSPR2 "
#endif

ساكن صحيح (*spr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr2_U, qspr2_L,
#وإلاإذا defined(DOUBLE)
  dspr2_U, dspr2_L,
#אחרת
  sspr2_U, sspr2_L,
#endif
};

#ifdef SMP
ساكن صحيح (*spr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, صحيح) = {
#ifdef XDOUBLE
  qspr2_thread_U, qspr2_thread_L,
#وإلاإذا defined(DOUBLE)
  dspr2_thread_U, dspr2_thread_L,
#אחרת
  sspr2_thread_U, sspr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  صحيح uplo;
  FLOAT *buffer;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  אם (uplo_arg  == 'U') uplo  = 0;
  אם (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  אם (incy == 0)          info =  7;
  אם (incx == 0)          info =  5;
  אם (n < 0)              info =  2;
  אם (uplo  < 0)          info =  1;

  אם (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#אחרת

فراغ CNAME(تعداد CBLAS_ORDER order,
	   تعداد CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a) {

  FLOAT *buffer;
  صحيح uplo;
  blasint info;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  אם (order == CblasColMajor) {
    אם (Uplo == CblasUpper)         uplo  = 0;
    אם (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    אם (incy == 0)          info =  7;
    אם (incx == 0)          info =  5;
    אם (n < 0)              info =  2;
    אם (uplo  < 0)          info =  1;
  }

  אם (order == CblasRowMajor) {
    אם (Uplo == CblasUpper)         uplo  = 1;
    אם (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    אם (incy == 0)          info =  7;
    אם (incx == 0)          info =  5;
    אם (n < 0)              info =  2;
    אם (uplo  < 0)          info =  1;
  }

  אם (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#endif

  אם (n == 0) החזר;

  אם (alpha == ZERO) החזר;

  אם (incx == 1 && incy == 1 && n < 50) {
    blasint i;
    אם (!uplo) {
      עבור (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, عدم, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, عدم, 0);
        a += i + 1;
      }
    } אחרת {
      עבור (i = 0; i < n; i++){
	AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, عدم, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, عدم, 0);
        a += n - i;
      }
    }
    החזר;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  אם (incx < 0 ) x -= (n - 1) * incx;
  אם (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  אם (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    החזר;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  אם (n < 50) nthreads = 1;
  אם (nthreads > 2 && n < 500) nthreads = 2;

  אם (nthreads == 1) {
#endif

  (spr2[uplo])(n, alpha, x, incx, y, incy, a, buffer);

#ifdef SMP
  } אחרת {

    (spr2_thread[uplo])(n, alpha, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  החזר;
}
