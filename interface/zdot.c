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

#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ RETURN_BY_STRUCT
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MYTYPE myxcomplex_t
#وگرنہ_اگر defined DOUBLE
#تعریف_کرو MYTYPE myzcomplex_t
#ورنہ
#تعریف_کرو MYTYPE myccomplex_t
#اگر_ختم
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

#اگر_تعریف_شدہ RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وگرنہ_اگر defined RETURN_BY_ڈھیر
خالی           NAME(OPENBLAS_COMPLEX_FLOAT *نتیجہ, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ورنہ
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#اگر_ختم

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#اگر_غیر_تعریف_شدہ RETURN_BY_ڈھیر
  OPENBLAS_COMPLEX_FLOAT ret;
#اگر_ختم
#اگر_تعریف_شدہ RETURN_BY_STRUCT
  MYTYPE  myret;
#اگر_ختم

#اگر_غیر_تعریف_شدہ RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#اگر_ختم

  PRINT_DEBUG_NAME;

  اگر (n <= 0) {
#اگر_تعریف_شدہ RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    واپس myret;
#وگرنہ_اگر defined RETURN_BY_ڈھیر
    *نتیجہ = zero;
    واپس;
#ورنہ
    واپس zero;
#اگر_ختم
  }

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0) x -= (n - 1) * incx * 2;
  اگر (incy < 0) y -= (n - 1) * incy * 2;

#اگر_تعریف_شدہ RETURN_BY_STRUCT

#اگر_غیر_تعریف_شدہ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ورنہ
  ret = DOTC_K(n, x, incx, y, incy);
#اگر_ختم

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

  واپس myret;

#وگرنہ_اگر defined RETURN_BY_ڈھیر

#اگر_غیر_تعریف_شدہ CONJ
  *نتیجہ = DOTU_K(n, x, incx, y, incy);
#ورنہ
  *نتیجہ = DOTC_K(n, x, incx, y, incy);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

#ورنہ

#اگر_غیر_تعریف_شدہ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ورنہ
  ret = DOTC_K(n, x, incx, y, incy);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

  واپس ret;

#اگر_ختم

}

#ورنہ

#اگر_تعریف_شدہ FORCE_USE_ڈھیر
خالی           CNAME(blasint n, خالی *vx, blasint incx, خالی *vy, blasint incy, خالی* vresult){
OPENBLAS_COMPLEX_FLOAT *نتیجہ= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#ورنہ
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, خالی *vx, blasint incx, خالی *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#اگر_ختم
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  اگر (n <= 0) {
#اگر_تعریف_شدہ FORCE_USE_ڈھیر
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *نتیجہ = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    واپس;
#ورنہ
    واپس zero;
#اگر_ختم
  }

  اگر (incx < 0) x -= (n - 1) * incx * 2;
  اگر (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

#اگر_تعریف_شدہ FORCE_USE_ڈھیر

#اگر_غیر_تعریف_شدہ CONJ
  *نتیجہ = DOTU_K(n, x, incx, y, incy);
#ورنہ
  *نتیجہ = DOTC_K(n, x, incx, y, incy);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

#ورنہ

#اگر_غیر_تعریف_شدہ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ورنہ
  ret = DOTC_K(n, x, incx, y, incy);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

  واپس ret;

#اگر_ختم

}

#اگر_ختم
