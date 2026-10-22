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
#تعریف_کړه ERROR_NAME "XHPR  "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZHPR  "
#وإلا
#تعریف_کړه ERROR_NAME "CHPR  "
#د_که_پای

ساكن صحيح (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#که_تعریف_شوی XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#وإلاإذا defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#وإلا
  chpr_U, chpr_L, chpr_V, chpr_M,
#د_که_پای
};

#که_تعریف_شوی SMP
ساكن صحيح (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, صحيح) = {
#که_تعریف_شوی XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#وإلاإذا defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#وإلا
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#د_که_پای
};
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

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

    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    إذا (Uplo == Cblas_پاسنی)         uplo  = 3;
    إذا (Uplo == Cblas_ښکتنی)         uplo  = 2;

    info = -1;

    إذا (incx == 0)          info =  5;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#د_که_پای

  إذا (n == 0) أرجع;

  إذا (alpha == ZERO) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0 ) x -= (n - 1) * incx * 2;

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

  (hpr[uplo])(n, alpha, x, incx, a, بفر);

#که_تعریف_شوی SMP

  } وإلا {

    (hpr_thread[uplo])(n, alpha, x, incx, a, بفر, nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(4, n * n / 2 + n, n * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
