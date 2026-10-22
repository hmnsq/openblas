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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XHPR2 "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZHPR2 "
#وإلا
#تعریف_کړه ERROR_NAME "CHPR2 "
#د_که_پای

ساكن صحيح (*hpr2[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#که_تعریف_شوی XDOUBLE
  xhpr2_U, xhpr2_L, xhpr2_V, xhpr2_M,
#وإلاإذا defined(DOUBLE)
  zhpr2_U, zhpr2_L, zhpr2_V, zhpr2_M,
#وإلا
  chpr2_U, chpr2_L, chpr2_V, chpr2_M,
#د_که_پای
};

#که_تعریف_شوی SMP
ساكن صحيح (*hpr2_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, صحيح) = {
#که_تعریف_شوی XDOUBLE
  xhpr2_thread_U, xhpr2_thread_L, xhpr2_thread_V, xhpr2_thread_M,
#وإلاإذا defined(DOUBLE)
  zhpr2_thread_U, zhpr2_thread_L, zhpr2_thread_V, zhpr2_thread_M,
#وإلا
  chpr2_thread_U, chpr2_thread_L, chpr2_thread_V, chpr2_thread_M,
#د_که_پای
};
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  صحيح uplo;
  FLOAT *بفر;
#که_تعریف_شوی SMP
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (incy == 0)          info =  7;
  إذا (incx == 0)          info =  5;
  إذا (n < 0)              info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   تعداد CBLAS_UPLO Uplo,
	   blasint n,
	   فراغ *VALPHA,
	   فراغ  *vx, blasint incx,
	   فراغ  *vy, blasint incy,
	   فراغ  *va) {

  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;
  FLOAT* a = (FLOAT*) va;

  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  FLOAT *بفر;
  صحيح uplo;
  blasint info;
#که_تعریف_شوی SMP
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_پاسنی)         uplo  = 0;
    إذا (Uplo == Cblas_ښکتنی)         uplo  = 1;

    info = -1;

    إذا (incy == 0)          info =  7;
    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    إذا (Uplo == Cblas_پاسنی)         uplo  = 3;
    إذا (Uplo == Cblas_ښکتنی)         uplo  = 2;

    info = -1;

    إذا (incx == 0)          info =  7;
    إذا (incy == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#د_که_پای

  إذا (n == 0) أرجع;

  إذا ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;


  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;
  إذا (incy < 0 ) y -= (n - 1) * incy * 2;

  بفر = (FLOAT *)blas_حافظه_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#که_تعریف_شوی SMP
  إذا (n < 100)
	nthreads = 1;
  وإلا
    nthreads = عدد_cpu_avail(2);

  إذا (nthreads == 1) {
#د_که_پای

    (hpr2[uplo])(n, alpha_r, alpha_i, x, incx, y, incy, a, بفر);

#که_تعریف_شوی SMP
  } وإلا {

    (hpr2_thread[uplo])(n, ALPHA, x, incx, y, incy, a, بفر, nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(4, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
