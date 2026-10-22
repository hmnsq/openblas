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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی RETURN_BY_STRUCT
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MYTYPE myxcomplex_t
#وإلاإذا defined DOUBLE
#تعریف_کړه MYTYPE myzcomplex_t
#وإلا
#تعریف_کړه MYTYPE myccomplex_t
#د_که_پای
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

#که_تعریف_شوی RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلاإذا defined RETURN_BY_ډېری
فراغ           NAME(OPENBLAS_COMPLEX_FLOAT *پایله, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#وإلا
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#د_که_پای

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#که_نه_تعریف_شوی RETURN_BY_ډېری
  OPENBLAS_COMPLEX_FLOAT ret;
#د_که_پای
#که_تعریف_شوی RETURN_BY_STRUCT
  MYTYPE  myret;
#د_که_پای

#که_نه_تعریف_شوی RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#د_که_پای

  PRINT_DEBUG_NAME;

  إذا (n <= 0) {
#که_تعریف_شوی RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    أرجع myret;
#وإلاإذا defined RETURN_BY_ډېری
    *پایله = zero;
    أرجع;
#وإلا
    أرجع zero;
#د_که_پای
  }

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#که_تعریف_شوی RETURN_BY_STRUCT

#که_نه_تعریف_شوی CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#د_که_پای

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_پای_ټکی(4, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع myret;

#وإلاإذا defined RETURN_BY_ډېری

#که_نه_تعریف_شوی CONJ
  *پایله = DOTU_K(n, x, incx, y, incy);
#وإلا
  *پایله = DOTC_K(n, x, incx, y, incy);
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(4, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

#وإلا

#که_نه_تعریف_شوی CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(4, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع ret;

#د_که_پای

}

#وإلا

#که_تعریف_شوی FORCE_USE_ډېری
فراغ           CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy, فراغ* vresult){
OPENBLAS_COMPLEX_FLOAT *پایله= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#وإلا
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#د_که_پای
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) {
#که_تعریف_شوی FORCE_USE_ډېری
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *پایله = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    أرجع;
#وإلا
    أرجع zero;
#د_که_پای
  }

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

#که_تعریف_شوی FORCE_USE_ډېری

#که_نه_تعریف_شوی CONJ
  *پایله = DOTU_K(n, x, incx, y, incy);
#وإلا
  *پایله = DOTC_K(n, x, incx, y, incy);
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(4, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

#وإلا

#که_نه_تعریف_شوی CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#وإلا
  ret = DOTC_K(n, x, incx, y, incy);
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(4, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع ret;

#د_که_پای

}

#د_که_پای
