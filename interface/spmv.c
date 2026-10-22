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
#تعریف_کړه ERROR_NAME "QSPMV "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DSPMV "
#وإلا
#تعریف_کړه ERROR_NAME "SSPMV "
#د_که_پای

ساكن  صحيح (*spmv[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#که_تعریف_شوی XDOUBLE
  qspmv_U, qspmv_L,
#وإلاإذا defined(DOUBLE)
  dspmv_U, dspmv_L,
#وإلا
  sspmv_U, sspmv_L,
#د_که_پای
};

#که_تعریف_شوی SMPTEST
ساكن  صحيح (*spmv_thread[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#که_تعریف_شوی XDOUBLE
  qspmv_thread_U, qspmv_thread_L,
#وإلاإذا defined(DOUBLE)
  dspmv_thread_U, dspmv_thread_L,
#وإلا
  sspmv_thread_U, sspmv_thread_L,
#د_که_پای
};
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  صحيح uplo;
  FLOAT *بفر;
#که_تعریف_شوی SMPTEST
  صحيح nthreads;
#د_که_پای

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

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   تعداد CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *a,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *بفر;
  صحيح uplo;
  blasint info;
#که_تعریف_شوی SMPTEST
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Uplo == Cblas_پاسنی)         uplo  = 0;
    إذا (Uplo == Cblas_ښکتنی)         uplo  = 1;

    info = -1;

    إذا (incy == 0)          info =  9;
    إذا (incx == 0)          info =  6;
    إذا (n < 0)              info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    إذا (Uplo == Cblas_پاسنی)         uplo  = 1;
    إذا (Uplo == Cblas_ښکتنی)         uplo  = 0;

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

#د_که_پای

  إذا (n == 0) أرجع;

  إذا (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);

  إذا (alpha == ZERO) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0 ) x -= (n - 1) * incx;
  إذا (incy < 0 ) y -= (n - 1) * incy;

  بفر = (FLOAT *)blas_حافظه_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#که_تعریف_شوی SMPTEST
  nthreads = عدد_cpu_avail(2);

  إذا (nthreads == 1) {
#د_که_پای

  (spmv[uplo])(n, alpha, a, x, incx, y, incy, بفر);

#که_تعریف_شوی SMPTEST
  } وإلا {

    (spmv_thread[uplo])(n, alpha, a, x, incx, y, incy, بفر, nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(1, n * n / 2 + n, n * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
