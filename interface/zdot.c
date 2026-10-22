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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef RETURN_BY_STRUCT
#ifdef XDOUBLE
#define MYTYPE myxcomplex_t
#ನೋಚೇದ್ಯದಿ defined DOUBLE
#define MYTYPE myzcomplex_t
#ಅನ್ಯಥಾ
#define MYTYPE myccomplex_t
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

#ifdef RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ನೋಚೇದ್ಯದಿ defined RETURN_BY_ರಾಶಿ
ಶೂನ್ಯ           NAME(OPENBLAS_COMPLEX_FLOAT *result, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ಅನ್ಯಥಾ
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#endif

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#ifndef RETURN_BY_ರಾಶಿ
  OPENBLAS_COMPLEX_FLOAT ret;
#endif
#ifdef RETURN_BY_STRUCT
  MYTYPE  myret;
#endif

#ifndef RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#endif

  PRINT_DEBUG_NAME;

  ಯದಿ (n <= 0) {
#ifdef RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    ಪ್ರತ್ಯಾಗಮ myret;
#ನೋಚೇದ್ಯದಿ defined RETURN_BY_ರಾಶಿ
    *result = zero;
    ಪ್ರತ್ಯಾಗಮ;
#ಅನ್ಯಥಾ
    ಪ್ರತ್ಯಾಗಮ zero;
#endif
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (n - 1) * incx * 2;
  ಯದಿ (incy < 0) y -= (n - 1) * incy * 2;

#ifdef RETURN_BY_STRUCT

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಅನ್ಯಥಾ
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ myret;

#ನೋಚೇದ್ಯದಿ defined RETURN_BY_ರಾಶಿ

#ifndef CONJ
  *result = DOTU_K(n, x, incx, y, incy);
#ಅನ್ಯಥಾ
  *result = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#ಅನ್ಯಥಾ

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಅನ್ಯಥಾ
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;

#endif

}

#ಅನ್ಯಥಾ

#ifdef FORCE_USE_ರಾಶಿ
ಶೂನ್ಯ           CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx, ಶೂನ್ಯ *vy, blasint incy, ಶೂನ್ಯ* vresult){
OPENBLAS_COMPLEX_FLOAT *result= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#ಅನ್ಯಥಾ
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx, ಶೂನ್ಯ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#endif
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  ಯದಿ (n <= 0) {
#ifdef FORCE_USE_ರಾಶಿ
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *result = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    ಪ್ರತ್ಯಾಗಮ;
#ಅನ್ಯಥಾ
    ಪ್ರತ್ಯಾಗಮ zero;
#endif
  }

  ಯದಿ (incx < 0) x -= (n - 1) * incx * 2;
  ಯದಿ (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifdef FORCE_USE_ರಾಶಿ

#ifndef CONJ
  *result = DOTU_K(n, x, incx, y, incy);
#ಅನ್ಯಥಾ
  *result = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#ಅನ್ಯಥಾ

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಅನ್ಯಥಾ
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;

#endif

}

#endif
