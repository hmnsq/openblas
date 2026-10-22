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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् RETURN_BY_STRUCT
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MYTYPE myxcomplex_t
#नोचेद्यदि defined DOUBLE
#परिभाषय MYTYPE myzcomplex_t
#अन्यथा
#परिभाषय MYTYPE myccomplex_t
#यदि_अन्तः
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

#यदि_परिभाषितम् RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#नोचेद्यदि defined RETURN_BY_राशिः
शून्य           NAME(OPENBLAS_COMPLEX_FLOAT *परिणाम, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#अन्यथा
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#यदि_अन्तः

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#यदि_अपरिभाषितम् RETURN_BY_राशिः
  OPENBLAS_COMPLEX_FLOAT ret;
#यदि_अन्तः
#यदि_परिभाषितम् RETURN_BY_STRUCT
  MYTYPE  myret;
#यदि_अन्तः

#यदि_अपरिभाषितम् RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  यदि (n <= 0) {
#यदि_परिभाषितम् RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    प्रत्यागम myret;
#नोचेद्यदि defined RETURN_BY_राशिः
    *परिणाम = zero;
    प्रत्यागम;
#अन्यथा
    प्रत्यागम zero;
#यदि_अन्तः
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (n - 1) * incx * 2;
  यदि (incy < 0) y -= (n - 1) * incy * 2;

#यदि_परिभाषितम् RETURN_BY_STRUCT

#यदि_अपरिभाषितम् CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#यदि_अन्तः

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

  प्रत्यागम myret;

#नोचेद्यदि defined RETURN_BY_राशिः

#यदि_अपरिभाषितम् CONJ
  *परिणाम = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  *परिणाम = DOTC_K(n, x, incx, y, incy);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

#अन्यथा

#यदि_अपरिभाषितम् CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

  प्रत्यागम ret;

#यदि_अन्तः

}

#अन्यथा

#यदि_परिभाषितम् FORCE_USE_राशिः
शून्य           CNAME(blasint n, शून्य *vx, blasint incx, शून्य *vy, blasint incy, शून्य* vresult){
OPENBLAS_COMPLEX_FLOAT *परिणाम= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#अन्यथा
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, शून्य *vx, blasint incx, शून्य *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#यदि_अन्तः
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  यदि (n <= 0) {
#यदि_परिभाषितम् FORCE_USE_राशिः
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *परिणाम = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    प्रत्यागम;
#अन्यथा
    प्रत्यागम zero;
#यदि_अन्तः
  }

  यदि (incx < 0) x -= (n - 1) * incx * 2;
  यदि (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#यदि_परिभाषितम् FORCE_USE_राशिः

#यदि_अपरिभाषितम् CONJ
  *परिणाम = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  *परिणाम = DOTC_K(n, x, incx, y, incy);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

#अन्यथा

#यदि_अपरिभाषितम् CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

  प्रत्यागम ret;

#यदि_अन्तः

}

#यदि_अन्तः
