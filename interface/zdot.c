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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MYTYPE myxcomplex_t
#ਜਾਂਜੇ defined DOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MYTYPE myzcomplex_t
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MYTYPE myccomplex_t
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ਜਾਂਜੇ defined RETURN_BY_ਢੇਰ
ਖਾਲੀ           NAME(OPENBLAS_COMPLEX_FLOAT *ਨਤੀਜਾ, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ਨਹੀਂਤਾਂ
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ਜੇ_ਅੰਤ

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ RETURN_BY_ਢੇਰ
  OPENBLAS_COMPLEX_FLOAT ret;
#ਜੇ_ਅੰਤ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT
  MYTYPE  myret;
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  ਜੇ (n <= 0) {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    ਵਾਪਸ myret;
#ਜਾਂਜੇ defined RETURN_BY_ਢੇਰ
    *ਨਤੀਜਾ = zero;
    ਵਾਪਸ;
#ਨਹੀਂਤਾਂ
    ਵਾਪਸ zero;
#ਜੇ_ਅੰਤ
  }

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (n - 1) * incx * 2;
  ਜੇ (incy < 0) y -= (n - 1) * incy * 2;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ RETURN_BY_STRUCT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ਨਹੀਂਤਾਂ
  ret = DOTC_K(n, x, incx, y, incy);
#ਜੇ_ਅੰਤ

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_ਸਿਰਾ(4, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ myret;

#ਜਾਂਜੇ defined RETURN_BY_ਢੇਰ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CONJ
  *ਨਤੀਜਾ = DOTU_K(n, x, incx, y, incy);
#ਨਹੀਂਤਾਂ
  *ਨਤੀਜਾ = DOTC_K(n, x, incx, y, incy);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(4, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

#ਨਹੀਂਤਾਂ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ਨਹੀਂਤਾਂ
  ret = DOTC_K(n, x, incx, y, incy);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(4, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;

#ਜੇ_ਅੰਤ

}

#ਨਹੀਂਤਾਂ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FORCE_USE_ਢੇਰ
ਖਾਲੀ           CNAME(blasint n, ਖਾਲੀ *vx, blasint incx, ਖਾਲੀ *vy, blasint incy, ਖਾਲੀ* vresult){
OPENBLAS_COMPLEX_FLOAT *ਨਤੀਜਾ= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#ਨਹੀਂਤਾਂ
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, ਖਾਲੀ *vx, blasint incx, ਖਾਲੀ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#ਜੇ_ਅੰਤ
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  ਜੇ (n <= 0) {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FORCE_USE_ਢੇਰ
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *ਨਤੀਜਾ = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    ਵਾਪਸ;
#ਨਹੀਂਤਾਂ
    ਵਾਪਸ zero;
#ਜੇ_ਅੰਤ
  }

  ਜੇ (incx < 0) x -= (n - 1) * incx * 2;
  ਜੇ (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FORCE_USE_ਢੇਰ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CONJ
  *ਨਤੀਜਾ = DOTU_K(n, x, incx, y, incy);
#ਨਹੀਂਤਾਂ
  *ਨਤੀਜਾ = DOTC_K(n, x, incx, y, incy);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(4, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

#ਨਹੀਂਤਾਂ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ਨਹੀਂਤਾਂ
  ret = DOTC_K(n, x, incx, y, incy);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(4, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;

#ਜੇ_ਅੰਤ

}

#ਜੇ_ਅੰਤ
