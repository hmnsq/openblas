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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#لغو_تعریف MAX_K

#اگر_تعریف_شده USE_ABS

#إذا defined(DOUBLE)
#تعریف_کن ABS fabs
#وإلا
#تعریف_کن ABS fabsf
#پایان_اگر

#اگر_تعریف_نشده USE_MIN

/* ABS & MAX */
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	QAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	DAMAX_K
#وإلا
#تعریف_کن MAX_K	SAMAX_K
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	XAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	ZAMAX_K
#وإلا
#تعریف_کن MAX_K	CAMAX_K
#پایان_اگر
#پایان_اگر

#وإلا

/* ABS & MIN */
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	QAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	DAMIN_K
#وإلا
#تعریف_کن MAX_K	SAMIN_K
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	XAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	ZAMIN_K
#وإلا
#تعریف_کن MAX_K	CAMIN_K
#پایان_اگر
#پایان_اگر

#پایان_اگر

#وإلا

#تعریف_کن ABS

#اگر_تعریف_نشده USE_MIN

/* MAX */
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	QMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	DMAX_K
#وإلا
#تعریف_کن MAX_K	SMAX_K
#پایان_اگر

#وإلا

/* MIN */
#اگر_تعریف_شده XDOUBLE
#تعریف_کن MAX_K	QMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کن MAX_K	DMIN_K
#وإلا
#تعریف_کن MAX_K	SMIN_K
#پایان_اگر

#پایان_اگر

#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع 0;

#اگر_تعریف_نشده COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#پایان_اگر

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_انتها(COMPSIZE, n, 0);

  IDEBUG_انتها;

  أرجع ret;
}

#وإلا

#اگر_تعریف_شده COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#پایان_اگر
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0;

#اگر_تعریف_نشده COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#پایان_اگر

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_انتها(COMPSIZE, n, 0);

  IDEBUG_انتها;

  أرجع ret;
}

#پایان_اگر
