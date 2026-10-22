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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_تعریف_شدہ RETURN_BY_STRUCT
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MYTYPE myxcomplex_t
#وإلاإذا defined DOUBLE
#تعریف_کرو MYTYPE myzcomplex_t
#وإلا
#تعریف_کرو MYTYPE myccomplex_t
#جے_انت
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

#جے_تعریف_شدہ RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلاإذا defined RETURN_BY_ڈھیر
فراغ           NAME(OPENBLAS_COMPLEX_FLOAT *نتیجہ, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلا
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#جے_انت

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#جے_بنا_تعریف RETURN_BY_ڈھیر
  OPENBLAS_COMPLEX_FLOAT ret;
#جے_انت
#جے_تعریف_شدہ RETURN_BY_STRUCT
  MYTYPE  myret;
#جے_انت

#جے_بنا_تعریف RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#جے_انت

  PRINT_DEBUG_NAME;

  إذا (n <= 0) {
#جے_تعریف_شدہ RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    أرجع myret;
#وإلاإذا defined RETURN_BY_ڈھیر
    *نتیجہ = zero;
    أرجع;
#وإلا
    أرجع zero;
#جے_انت
  }

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#جے_تعریف_شدہ RETURN_BY_STRUCT

#جے_بنا_تعریف CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جے_انت

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

  أرجع myret;

#وإلاإذا defined RETURN_BY_ڈھیر

#جے_بنا_تعریف CONJ
  *نتیجہ = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتیجہ = DOTC_K(n, x, incx, y, incy);
#جے_انت

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

#وإلا

#جے_بنا_تعریف CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جے_انت

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

  أرجع ret;

#جے_انت

}

#وإلا

#جے_تعریف_شدہ FORCE_USE_ڈھیر
فراغ           CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy, فراغ* vresult){
OPENBLAS_COMPLEX_FLOAT *نتیجہ= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#وإلا
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#جے_انت
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) {
#جے_تعریف_شدہ FORCE_USE_ڈھیر
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *نتیجہ = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    أرجع;
#وإلا
    أرجع zero;
#جے_انت
  }

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#جے_تعریف_شدہ FORCE_USE_ڈھیر

#جے_بنا_تعریف CONJ
  *نتیجہ = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتیجہ = DOTC_K(n, x, incx, y, incy);
#جے_انت

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

#وإلا

#جے_بنا_تعریف CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جے_انت

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

  أرجع ret;

#جے_انت

}

#جے_انت
