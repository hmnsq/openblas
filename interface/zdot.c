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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل RETURN_BY_STRUCT
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MYTYPE myxcomplex_t
#وإلاإذا defined DOUBLE
#وصف_ڪريو MYTYPE myzcomplex_t
#وإلا
#وصف_ڪريو MYTYPE myccomplex_t
#جيڪڏهن_پورو
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

#جيڪڏهن_وصف_ٿيل RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلاإذا defined RETURN_BY_ڍير
فراغ           NAME(OPENBLAS_COMPLEX_FLOAT *نتيجو, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلا
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#جيڪڏهن_پورو

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#جيڪڏهن_وصف_نه_ٿيل RETURN_BY_ڍير
  OPENBLAS_COMPLEX_FLOAT ret;
#جيڪڏهن_پورو
#جيڪڏهن_وصف_ٿيل RETURN_BY_STRUCT
  MYTYPE  myret;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#جيڪڏهن_پورو

  PRINT_DEBUG_NAME;

  إذا (n <= 0) {
#جيڪڏهن_وصف_ٿيل RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    أرجع myret;
#وإلاإذا defined RETURN_BY_ڍير
    *نتيجو = zero;
    أرجع;
#وإلا
    أرجع zero;
#جيڪڏهن_پورو
  }

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#جيڪڏهن_وصف_ٿيل RETURN_BY_STRUCT

#جيڪڏهن_وصف_نه_ٿيل CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جيڪڏهن_پورو

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

  أرجع myret;

#وإلاإذا defined RETURN_BY_ڍير

#جيڪڏهن_وصف_نه_ٿيل CONJ
  *نتيجو = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتيجو = DOTC_K(n, x, incx, y, incy);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

#وإلا

#جيڪڏهن_وصف_نه_ٿيل CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

  أرجع ret;

#جيڪڏهن_پورو

}

#وإلا

#جيڪڏهن_وصف_ٿيل FORCE_USE_ڍير
فراغ           CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy, فراغ* vresult){
OPENBLAS_COMPLEX_FLOAT *نتيجو= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#وإلا
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#جيڪڏهن_پورو
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) {
#جيڪڏهن_وصف_ٿيل FORCE_USE_ڍير
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *نتيجو = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    أرجع;
#وإلا
    أرجع zero;
#جيڪڏهن_پورو
  }

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

#جيڪڏهن_وصف_ٿيل FORCE_USE_ڍير

#جيڪڏهن_وصف_نه_ٿيل CONJ
  *نتيجو = DOTU_K(n, x, incx, y, incy);
#وإلا
  *نتيجو = DOTC_K(n, x, incx, y, incy);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

#وإلا

#جيڪڏهن_وصف_نه_ٿيل CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

  أرجع ret;

#جيڪڏهن_پورو

}

#جيڪڏهن_پورو
