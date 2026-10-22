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
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef RETURN_BY_STRUCT
#ifdef XDOUBLE
#define MYTYPE myxcomplex_t
#नोचेद्यदि defined DOUBLE
#define MYTYPE myzcomplex_t
#अञ्ञथा
#define MYTYPE myccomplex_t
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

#ifdef RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#नोचेद्यदि defined RETURN_BY_राशिः
शून्य           NAME(OPENBLAS_COMPLEX_FLOAT *result, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#अञ्ञथा
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#endif

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#ifndef RETURN_BY_राशिः
  OPENBLAS_COMPLEX_FLOAT ret;
#endif
#ifdef RETURN_BY_STRUCT
  MYTYPE  myret;
#endif

#ifndef RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#endif

  PRINT_DEBUG_NAME;

  सचे (n <= 0) {
#ifdef RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    पटिनिवत्त myret;
#नोचेद्यदि defined RETURN_BY_राशिः
    *result = zero;
    पटिनिवत्त;
#अञ्ञथा
    पटिनिवत्त zero;
#endif
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0) x -= (n - 1) * incx * 2;
  सचे (incy < 0) y -= (n - 1) * incy * 2;

#ifdef RETURN_BY_STRUCT

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अञ्ञथा
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  पटिनिवत्त myret;

#नोचेद्यदि defined RETURN_BY_राशिः

#ifndef CONJ
  *result = DOTU_K(n, x, incx, y, incy);
#अञ्ञथा
  *result = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#अञ्ञथा

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अञ्ञथा
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  पटिनिवत्त ret;

#endif

}

#अञ्ञथा

#ifdef FORCE_USE_राशिः
शून्य           CNAME(blasint n, शून्य *vx, blasint incx, शून्य *vy, blasint incy, शून्य* vresult){
OPENBLAS_COMPLEX_FLOAT *result= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#अञ्ञथा
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, शून्य *vx, blasint incx, शून्य *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#endif
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  सचे (n <= 0) {
#ifdef FORCE_USE_राशिः
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *result = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    पटिनिवत्त;
#अञ्ञथा
    पटिनिवत्त zero;
#endif
  }

  सचे (incx < 0) x -= (n - 1) * incx * 2;
  सचे (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifdef FORCE_USE_राशिः

#ifndef CONJ
  *result = DOTU_K(n, x, incx, y, incy);
#अञ्ञथा
  *result = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

#अञ्ञथा

#ifndef CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#अञ्ञथा
  ret = DOTC_K(n, x, incx, y, incy);
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  पटिनिवत्त ret;

#endif

}

#endif
