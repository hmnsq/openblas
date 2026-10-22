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
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_شده RETURN_BY_STRUCT
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MYTYPE myxcomplex_t
#وإلاإذا defined DOUBLE
#تعریف_کن MYTYPE myzcomplex_t
#وإلا
#تعریف_کن MYTYPE myccomplex_t
#پایان_اگر
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

#اگر_تعریف_شده RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلاإذا defined RETURN_BY_پشته
فراغ           NAME(OPENBLAS_COMPLEX_FLOAT *نتیجه, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلا
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#پایان_اگر

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#اگر_تعریف_نشده RETURN_BY_پشته
  OPENBLAS_COMPLEX_FLOAT ret;
#پایان_اگر
#اگر_تعریف_شده RETURN_BY_STRUCT
  MYTYPE  myret;
#پایان_اگر

#اگر_تعریف_نشده RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#پایان_اگر

  PRINT_DEBUG_NAME;

  إذا (n <= 0) {
#اگر_تعریف_شده RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    أرجع myret;
#وإلاإذا defined RETURN_BY_پشته
    *نتیجه = zero;
    أرجع;
#وإلا
    أرجع zero;
#پایان_اگر
  }

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#اگر_تعریف_شده RETURN_BY_STRUCT

#اگر_تعریف_نشده CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#پایان_اگر

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_انتها(4, 2 * n, 2 * n);

  IDEBUG_انتها;

  أرجع myret;

#وإلاإذا defined RETURN_BY_پشته

#اگر_تعریف_نشده CONJ
  *نتیجه = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتیجه = DOTC_K(n, x, incx, y, incy);
#پایان_اگر

  FUNCTION_PROFILE_انتها(4, 2 * n, 2 * n);

  IDEBUG_انتها;

#وإلا

#اگر_تعریف_نشده CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#پایان_اگر

  FUNCTION_PROFILE_انتها(4, 2 * n, 2 * n);

  IDEBUG_انتها;

  أرجع ret;

#پایان_اگر

}

#وإلا

#اگر_تعریف_شده FORCE_USE_پشته
فراغ           CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy, فراغ* vresult){
OPENBLAS_COMPLEX_FLOAT *نتیجه= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#وإلا
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#پایان_اگر
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) {
#اگر_تعریف_شده FORCE_USE_پشته
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *نتیجه = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    أرجع;
#وإلا
    أرجع zero;
#پایان_اگر
  }

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#اگر_تعریف_شده FORCE_USE_پشته

#اگر_تعریف_نشده CONJ
  *نتیجه = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتیجه = DOTC_K(n, x, incx, y, incy);
#پایان_اگر

  FUNCTION_PROFILE_انتها(4, 2 * n, 2 * n);

  IDEBUG_انتها;

#وإلا

#اگر_تعریف_نشده CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#پایان_اگر

  FUNCTION_PROFILE_انتها(4, 2 * n, 2 * n);

  IDEBUG_انتها;

  أرجع ret;

#پایان_اگر

}

#پایان_اگر
