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

#تعریف_منسوخ MAX_K

#اگر_تعریف_شدہ USE_ABS

#اگر defined(DOUBLE)
#تعریف_کرو ABS fabs
#ورنہ
#تعریف_کرو ABS fabsf
#اگر_ختم

#اگر_غیر_تعریف_شدہ USE_MIN

/* ABS & MAX */
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QAMAX_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	DAMAX_K
#ورنہ
#تعریف_کرو MAX_K	SAMAX_K
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	XAMAX_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	ZAMAX_K
#ورنہ
#تعریف_کرو MAX_K	CAMAX_K
#اگر_ختم
#اگر_ختم

#ورنہ

/* ABS & MIN */
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QAMIN_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	DAMIN_K
#ورنہ
#تعریف_کرو MAX_K	SAMIN_K
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	XAMIN_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	ZAMIN_K
#ورنہ
#تعریف_کرو MAX_K	CAMIN_K
#اگر_ختم
#اگر_ختم

#اگر_ختم

#ورنہ

#تعریف_کرو ABS

#اگر_غیر_تعریف_شدہ USE_MIN

/* MAX */
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QMAX_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	DMAX_K
#ورنہ
#تعریف_کرو MAX_K	SMAX_K
#اگر_ختم

#ورنہ

/* MIN */
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو MAX_K	QMIN_K
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو MAX_K	DMIN_K
#ورنہ
#تعریف_کرو MAX_K	SMIN_K
#اگر_ختم

#اگر_ختم

#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  اگر (n <= 0) واپس 0;

#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (incx == 0) واپس (ABS(*x));
#ورنہ
  اگر (incx == 0) واپس (ABS(*x) + ABS(*(x+1)));
#اگر_ختم

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_انتہا(COMPSIZE, n, 0);

  IDEBUG_انتہا;

  واپس ret;
}

#ورنہ

#اگر_تعریف_شدہ COMPLEX
FLOAT CNAME(blasint n, خالی *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ورنہ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#اگر_ختم
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  اگر (n <= 0) واپس 0;

#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (incx == 0) واپس (ABS(*x));
#ورنہ
  اگر (incx == 0) واپس (ABS(*x) + ABS(*(x+1)));
#اگر_ختم

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_انتہا(COMPSIZE, n, 0);

  IDEBUG_انتہا;

  واپس ret;
}

#اگر_ختم
