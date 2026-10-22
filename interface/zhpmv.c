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
#وصف_ڪريو ERROR_NAME "XHPMV "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHPMV "
#وإلا
#وصف_ڪريو ERROR_NAME "CHPMV "
#جيڪڏهن_پورو

ساكن  صحيح (*hpmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  xhpmv_U, xhpmv_L, xhpmv_V, xhpmv_M,
#وإلاإذا defined(DOUBLE)
  zhpmv_U, zhpmv_L, zhpmv_V, zhpmv_M,
#وإلا
  chpmv_U, chpmv_L, chpmv_V, chpmv_M,
#جيڪڏهن_پورو
};

#جيڪڏهن_وصف_ٿيل SMP
ساكن  صحيح (*hpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  xhpmv_thread_U, xhpmv_thread_L, xhpmv_thread_V, xhpmv_thread_M,
#وإلاإذا defined(DOUBLE)
  zhpmv_thread_U, zhpmv_thread_L, zhpmv_thread_V, zhpmv_thread_M,
#وإلا
  chpmv_thread_U, chpmv_thread_L, chpmv_thread_V, chpmv_thread_M,
#جيڪڏهن_پورو
};
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  صحيح uplo;
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

  إذا (incy == 0)          info =  9;
  إذا (incx == 0)          info =  6;
  إذا (n < 0)              info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب,
	   تعداد CBLAS_UPLO Uplo,
	   blasint n,
	   فراغ *VALPHA,
	   فراغ  *va,
	   فراغ  *vx, blasint incx,
	   فراغ *VBETA,
	   فراغ  *vy, blasint incy){

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  FLOAT *بفر;
  صحيح uplo;
  blasint info;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (Uplo == Cblas_مٿانهون)         uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون)         uplo  = 1;

    info = -1;

    إذا (incy == 0)          info =  9;
    إذا (incx == 0)          info =  6;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (Uplo == Cblas_مٿانهون)         uplo  = 3;
    إذا (Uplo == Cblas_هيٺاهون)         uplo  = 2;

    info = -1;

    إذا (incy == 0)          info =  9;
    إذا (incx == 0)          info =  6;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_پورو

  إذا (n == 0) أرجع;

  إذا ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, y, blasabs(incy), عدم, 0, عدم, 0);

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
  nthreads = انگ_cpu_avail(2);

  إذا (n < 50) nthreads = 1;
  إذا (nthreads > 2 && n < 500) nthreads = 2;

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

  (hpmv[uplo])(n, alpha_r, alpha_i, a, x, incx, y, incy, بفر);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

    (hpmv_thread[uplo])(n, ALPHA, a, x, incx, y, incy, بفر, nthreads);

  }
#جيڪڏهن_پورو

  blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(4, n * n / 2 + n, n * n);

  IDEBUG_پڇاڙي;

  أرجع;
}
