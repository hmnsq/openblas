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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XHER2 "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHER2 "
#وإلا
#وصف_ڪريو ERROR_NAME "CHER2 "
#جيڪڏهن_پورو

ساكن صحيح (*her2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  xher2_U, xher2_L, xher2_V, xher2_M,
#وإلاإذا defined(DOUBLE)
  zher2_U, zher2_L, zher2_V, zher2_M,
#وإلا
  cher2_U, cher2_L, cher2_V, cher2_M,
#جيڪڏهن_پورو
};

#جيڪڏهن_وصف_ٿيل SMP
ساكن صحيح (*her2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  xher2_thread_U, xher2_thread_L, xher2_thread_V, xher2_thread_M,
#وإلاإذا defined(DOUBLE)
  zher2_thread_U, zher2_thread_L, zher2_thread_V, zher2_thread_M,
#وإلا
  cher2_thread_U, cher2_thread_L, cher2_thread_V, cher2_thread_M,
#جيڪڏهن_پورو
};
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  صحيح  uplo;
  FLOAT *بفر;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (lda  < MAX(1, n))   info =  9;
  إذا (incy == 0)          info =  7;
  إذا (incx == 0)          info =  5;
  إذا (n < 0)              info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_UPLO Uplo, blasint n, فراغ *VALPHA, فراغ *vx, blasint incx, فراغ *vy, blasint incy, فراغ *va, blasint lda) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT *بفر;
  صحيح uplo;
  blasint info;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {

    إذا (Uplo == Cblas_مٿانهون) uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 1;

    info = -1;

    إذا (lda  < MAX(1, n))   info =  9;
    إذا (incy == 0)          info =  7;
    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {

    إذا (Uplo == Cblas_مٿانهون) uplo  = 3;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 2;

    info = -1;

    إذا (lda  < MAX(1, n))   info =  9;
    إذا (incx == 0)          info =  7;
    إذا (incy == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_پورو

  إذا (n == 0) أرجع;

  إذا ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;
  إذا (incy < 0 ) y -= (n - 1) * incy * 2;

  بفر = (FLOAT *)blas_يادگيري_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_وصف_ٿيل SMP
  إذا (n < 100)
	nthreads = 1;
  وإلا
    nthreads = انگ_cpu_avail(2);

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

    (her2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, lda, بفر);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

    (her2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, lda, بفر, nthreads);

  }
#جيڪڏهن_پورو

  blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_پڇاڙي;

  أرجع;
}
