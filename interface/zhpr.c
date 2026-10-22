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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XHPR  "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZHPR  "
#وإلا
#تعریف_کن ERROR_NAME "CHPR  "
#پایان_اگر

ساكن صحيح (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#اگر_تعریف_شده XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#وإلاإذا defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#وإلا
  chpr_U, chpr_L, chpr_V, chpr_M,
#پایان_اگر
};

#اگر_تعریف_شده SMP
ساكن صحيح (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, صحيح) = {
#اگر_تعریف_شده XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#وإلاإذا defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#وإلا
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#پایان_اگر
};
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

  blasint info;
  صحيح uplo;
  FLOAT *میانگیر;
#اگر_تعریف_شده SMP
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

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
	   FLOAT alpha,
	   فراغ  *vx, blasint incx,
	   فراغ  *va) {

  FLOAT* x = (FLOAT*) vx;
  FLOAT* a = (FLOAT*) va;

  FLOAT *میانگیر;
  صحيح uplo;
  blasint info;
#اگر_تعریف_شده SMP
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_فوقانی)         uplo  = 0;
    إذا (Uplo == Cblas_پایینی)         uplo  = 1;

    info = -1;

    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {
    إذا (Uplo == Cblas_فوقانی)         uplo  = 3;
    إذا (Uplo == Cblas_پایینی)         uplo  = 2;

    info = -1;

    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#پایان_اگر

  إذا (n == 0) أرجع;

  إذا (alpha == ZERO) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;

  میانگیر = (FLOAT *)blas_حافظه_alloc(1);
  إذا (!میانگیر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#اگر_تعریف_شده SMP
  إذا (n < 100)
    nthreads = 1;
  وإلا
    nthreads = رقم_cpu_avail(2);

  إذا (nthreads == 1) {
#پایان_اگر

  (hpr[uplo])(n, alpha, x, incx, a, میانگیر);

#اگر_تعریف_شده SMP

  } وإلا {

    (hpr_thread[uplo])(n, alpha, x, incx, a, میانگیر, nthreads);

  }
#پایان_اگر

  blas_حافظه_free(میانگیر);

  FUNCTION_PROFILE_انتها(4, n * n / 2 + n, n * n);

  IDEBUG_انتها;

  أرجع;
}
