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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHER  "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHER  "
#وإلا
#تعریف_کرو ERROR_NAME "CHER  "
#جے_انت

ساكن صحيح (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#جے_تعریف_شدہ XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#وإلاإذا defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#وإلا
  cher_U, cher_L, cher_V, cher_M,
#جے_انت
};

#جے_تعریف_شدہ SMP
ساكن صحيح (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#جے_تعریف_شدہ XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#وإلاإذا defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#وإلا
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#جے_انت
};
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  صحيح uplo;
  FLOAT *بفر;
#جے_تعریف_شدہ SMP
  صحيح nthreads;
#جے_انت

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (lda  < MAX(1, n))   info =  7;
  إذا (incx == 0)          info =  5;
  إذا (n < 0)              info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo, blasint n, FLOAT alpha, فراغ *vx, blasint incx, فراغ *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *بفر;
  صحيح uplo;
  blasint info;
#جے_تعریف_شدہ SMP
  صحيح nthreads;
#جے_انت

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {

    إذا (Uplo == Cblas_اتلا) uplo  = 0;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 1;

    info = -1;

    إذا (lda  < MAX(1, n))   info =  7;
    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;

  }

  إذا (ترتیب == Cblas_صف_Major) {

    إذا (Uplo == Cblas_اتلا) uplo  = 3;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 2;

    info = -1;

    إذا (lda  < MAX(1, n))   info =  7;
    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جے_انت

  إذا (n == 0) أرجع;

  إذا (alpha == ZERO) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;

  بفر = (FLOAT *)blas_حافظہ_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جے_تعریف_شدہ SMP
  إذا (n < 100)
    nthreads = 1;
  وإلا
    nthreads = عدد_cpu_avail(2);

  إذا (nthreads == 1) {
#جے_انت

  (her[uplo])(n, alpha, x, incx, a, lda, بفر);

#جے_تعریف_شدہ SMP
  } وإلا {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, بفر, nthreads);

  }
#جے_انت

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_سرا(4, n * n / 2 + n, n * n);

  IDEBUG_سرا;

  أرجع;
}
