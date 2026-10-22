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

#تعریف_رد MAX_K

#جے_تعریف_شدہ USE_ABS

#إذا defined(DOUBLE)
#تعریف_کرو ABS fabs
#وإلا
#تعریف_کرو ABS fabsf
#جے_انت

#جے_بنا_تعریف USE_MIN

/* ABS & MAX */
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	DAMAX_K
#وإلا
#تعریف_کرو MAX_K	SAMAX_K
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	XAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	ZAMAX_K
#وإلا
#تعریف_کرو MAX_K	CAMAX_K
#جے_انت
#جے_انت

#وإلا

/* ABS & MIN */
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	DAMIN_K
#وإلا
#تعریف_کرو MAX_K	SAMIN_K
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	XAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	ZAMIN_K
#وإلا
#تعریف_کرو MAX_K	CAMIN_K
#جے_انت
#جے_انت

#جے_انت

#وإلا

#تعریف_کرو ABS

#جے_بنا_تعریف USE_MIN

/* MAX */
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	DMAX_K
#وإلا
#تعریف_کرو MAX_K	SMAX_K
#جے_انت

#وإلا

/* MIN */
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کرو MAX_K	DMIN_K
#وإلا
#تعریف_کرو MAX_K	SMIN_K
#جے_انت

#جے_انت

#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع 0;

#جے_بنا_تعریف COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#جے_انت

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_سرا(COMPSIZE, n, 0);

  IDEBUG_سرا;

  أرجع ret;
}

#وإلا

#جے_تعریف_شدہ COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#جے_انت
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0;

#جے_بنا_تعریف COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#جے_انت

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_سرا(COMPSIZE, n, 0);

  IDEBUG_سرا;

  أرجع ret;
}

#جے_انت
