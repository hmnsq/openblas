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

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  اگر (n <= 0) واپس 0.;

#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (n == 1)
#اگر_تعریف_شدہ DOUBLE
    واپس fabs(x[0]);
#ورنہ
    واپس fabsf(x[0]);
#اگر_ختم
#اگر_ختم

  اگر (incx == 0)
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ DOUBLE
  واپس (جذر((دوہرا)n)*fabs(x[0]));
#ورنہ
  واپس (جذر((اعشاری)n)*fabsf(x[0]));
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ DOUBLE
  {
  دوہرا fr=fabs(x[0]);
  دوہرا fi=fabs(x[1]);
  دوہرا fmin=MIN(fr,fi);
  دوہرا fmax=MAX(fr,fi);
  اگر (fmax==0.) واپس(fmax);
  اگر (fmax==fmin) واپس(جذر((دوہرا)n)*جذر(2.)*fmax);
  واپس (جذر((دوہرا)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ورنہ
  {
  اعشاری fr=fabs(x[0]);
  اعشاری fi=fabs(x[1]);
  اعشاری fmin=MIN(fr,fi);
  اعشاری fmax=MAX(fr,fi);
  اگر (fmax==0.) واپس(fmax);
  اگر (fmax==fmin) واپس(جذر((اعشاری)n)*جذر(2.)*fmax);
  واپس (جذر((اعشاری)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#اگر_ختم
#اگر_ختم

  اگر (incx < 0) 
#اگر_تعریف_شدہ COMPLEX    
    x -= (n - 1) * incx * 2;
#ورنہ
    x -= (n - 1) * incx;
#اگر_ختم
  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_انتہا(COMPSIZE, n, 2 * n);

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

  اگر (n <= 0) واپس 0.;

#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (n == 1)
#اگر_تعریف_شدہ DOUBLE
    واپس fabs(x[0]);
#ورنہ
    واپس fabsf(x[0]);
#اگر_ختم
#اگر_ختم

  اگر (incx == 0)
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ DOUBLE
  واپس (جذر((دوہرا)n)*fabs(x[0]));
#ورنہ
  واپس (جذر((اعشاری)n)*fabsf(x[0]));
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ DOUBLE
  {
  دوہرا fr=fabs(x[0]);
  دوہرا fi=fabs(x[1]);
  دوہرا fmin=MIN(fr,fi);
  دوہرا fmax=MAX(fr,fi);
  اگر (fmax==0.) واپس(fmax);
  اگر (fmax==fmin) واپس(جذر((دوہرا)n)*جذر(2.)*fmax);
  واپس (جذر((دوہرا)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ورنہ
  {
  اعشاری fr=fabs(x[0]);
  اعشاری fi=fabs(x[1]);
  اعشاری fmin=MIN(fr,fi);
  اعشاری fmax=MAX(fr,fi);
  اگر (fmax==0.) واپس(fmax);
  اگر (fmax==fmin) واپس(جذر((اعشاری)n)*جذر(2.)*fmax);
  واپس (جذر((اعشاری)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#اگر_ختم
#اگر_ختم

  اگر (incx < 0) 
#اگر_تعریف_شدہ COMPLEX    
    x -= (n - 1) * incx * 2;
#ورنہ
    x -= (n - 1) * incx;
#اگر_ختم
  
  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_انتہا(COMPSIZE, n, 2 * n);

  IDEBUG_انتہا;

  واپس ret;
}

#اگر_ختم
