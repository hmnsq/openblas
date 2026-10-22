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
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  اگر (n <= 0) بازگشت 0.;

#ifndef COMPLEX
  اگر (n == 1)
#ifdef DOUBLE
    بازگشت fabs(x[0]);
#وگرنه
    بازگشت fabsf(x[0]);
#endif
#endif

  اگر (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  بازگشت (sqrt((دوگانه)n)*fabs(x[0]));
#وگرنه
  بازگشت (sqrt((اعشاری)n)*fabsf(x[0]));
#endif
#وگرنه
#ifdef DOUBLE
  {
  دوگانه fr=fabs(x[0]);
  دوگانه fi=fabs(x[1]);
  دوگانه fmin=MIN(fr,fi);
  دوگانه fmax=MAX(fr,fi);
  اگر (fmax==0.) بازگشت(fmax);
  اگر (fmax==fmin) بازگشت(sqrt((دوگانه)n)*sqrt(2.)*fmax);
  بازگشت (sqrt((دوگانه)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وگرنه
  {
  اعشاری fr=fabs(x[0]);
  اعشاری fi=fabs(x[1]);
  اعشاری fmin=MIN(fr,fi);
  اعشاری fmax=MAX(fr,fi);
  اگر (fmax==0.) بازگشت(fmax);
  اگر (fmax==fmin) بازگشت(sqrt((اعشاری)n)*sqrt(2.)*fmax);
  بازگشت (sqrt((اعشاری)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  اگر (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#وگرنه
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  بازگشت ret;
}

#وگرنه

#ifdef COMPLEX
FLOAT CNAME(blasint n, تهی *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وگرنه
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  اگر (n <= 0) بازگشت 0.;

#ifndef COMPLEX
  اگر (n == 1)
#ifdef DOUBLE
    بازگشت fabs(x[0]);
#وگرنه
    بازگشت fabsf(x[0]);
#endif
#endif

  اگر (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  بازگشت (sqrt((دوگانه)n)*fabs(x[0]));
#وگرنه
  بازگشت (sqrt((اعشاری)n)*fabsf(x[0]));
#endif
#وگرنه
#ifdef DOUBLE
  {
  دوگانه fr=fabs(x[0]);
  دوگانه fi=fabs(x[1]);
  دوگانه fmin=MIN(fr,fi);
  دوگانه fmax=MAX(fr,fi);
  اگر (fmax==0.) بازگشت(fmax);
  اگر (fmax==fmin) بازگشت(sqrt((دوگانه)n)*sqrt(2.)*fmax);
  بازگشت (sqrt((دوگانه)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وگرنه
  {
  اعشاری fr=fabs(x[0]);
  اعشاری fi=fabs(x[1]);
  اعشاری fmin=MIN(fr,fi);
  اعشاری fmax=MAX(fr,fi);
  اگر (fmax==0.) بازگشت(fmax);
  اگر (fmax==fmin) بازگشت(sqrt((اعشاری)n)*sqrt(2.)*fmax);
  بازگشت (sqrt((اعشاری)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  اگر (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#وگرنه
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  بازگشت ret;
}

#endif
