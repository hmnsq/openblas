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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित RETURN_BY_STRUCT
#अगर_घोषित XDOUBLE
#घोषणा MYTYPE myxcomplex_t
#या_अगर defined DOUBLE
#घोषणा MYTYPE myzcomplex_t
#अन्यथा
#घोषणा MYTYPE myccomplex_t
#पूर्ण_अगर
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

#अगर_घोषित RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#या_अगर defined RETURN_BY_ढेर
व्योम           NAME(OPENBLAS_COMPLEX_FLOAT *परिणाम, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#अन्यथा
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#पूर्ण_अगर

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#अगर_अघोषित RETURN_BY_ढेर
  OPENBLAS_COMPLEX_FLOAT ret;
#पूर्ण_अगर
#अगर_घोषित RETURN_BY_STRUCT
  MYTYPE  myret;
#पूर्ण_अगर

#अगर_अघोषित RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  अगर (n <= 0) {
#अगर_घोषित RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    वापस myret;
#या_अगर defined RETURN_BY_ढेर
    *परिणाम = zero;
    वापस;
#अन्यथा
    वापस zero;
#पूर्ण_अगर
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (n - 1) * incx * 2;
  अगर (incy < 0) y -= (n - 1) * incy * 2;

#अगर_घोषित RETURN_BY_STRUCT

#अगर_अघोषित CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#पूर्ण_अगर

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  वापस myret;

#या_अगर defined RETURN_BY_ढेर

#अगर_अघोषित CONJ
  *परिणाम = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  *परिणाम = DOTC_K(n, x, incx, y, incy);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#अन्यथा

#अगर_अघोषित CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  वापस ret;

#पूर्ण_अगर

}

#अन्यथा

#अगर_घोषित FORCE_USE_ढेर
व्योम           CNAME(blasint n, व्योम *vx, blasint incx, व्योम *vy, blasint incy, व्योम* vresult){
OPENBLAS_COMPLEX_FLOAT *परिणाम= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#अन्यथा
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, व्योम *vx, blasint incx, व्योम *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#पूर्ण_अगर
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  अगर (n <= 0) {
#अगर_घोषित FORCE_USE_ढेर
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *परिणाम = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    वापस;
#अन्यथा
    वापस zero;
#पूर्ण_अगर
  }

  अगर (incx < 0) x -= (n - 1) * incx * 2;
  अगर (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#अगर_घोषित FORCE_USE_ढेर

#अगर_अघोषित CONJ
  *परिणाम = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  *परिणाम = DOTC_K(n, x, incx, y, incy);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#अन्यथा

#अगर_अघोषित CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अन्यथा
  ret = DOTC_K(n, x, incx, y, incy);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  वापस ret;

#पूर्ण_अगर

}

#पूर्ण_अगर
