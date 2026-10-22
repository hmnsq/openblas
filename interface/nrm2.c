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

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ਜੇ (n <= 0) ਵਾਪਸ 0.;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (n == 1)
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
    ਵਾਪਸ fabs(x[0]);
#ਨਹੀਂਤਾਂ
    ਵਾਪਸ fabsf(x[0]);
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  ਜੇ (incx == 0)
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n)*fabs(x[0]));
#ਨਹੀਂਤਾਂ
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n)*fabsf(x[0]));
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  {
  ਦੁੱਗਣਾ fr=fabs(x[0]);
  ਦੁੱਗਣਾ fi=fabs(x[1]);
  ਦੁੱਗਣਾ fmin=MIN(fr,fi);
  ਦੁੱਗਣਾ fmax=MAX(fr,fi);
  ਜੇ (fmax==0.) ਵਾਪਸ(fmax);
  ਜੇ (fmax==fmin) ਵਾਪਸ(ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n)*ਵਰਗਮੂਲ(2.)*fmax);
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n) * fmax * ਵਰਗਮੂਲ (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ਨਹੀਂਤਾਂ
  {
  ਦਸ਼ਮਲਵ fr=fabs(x[0]);
  ਦਸ਼ਮਲਵ fi=fabs(x[1]);
  ਦਸ਼ਮਲਵ fmin=MIN(fr,fi);
  ਦਸ਼ਮਲਵ fmax=MAX(fr,fi);
  ਜੇ (fmax==0.) ਵਾਪਸ(fmax);
  ਜੇ (fmax==fmin) ਵਾਪਸ(ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n)*ਵਰਗਮੂਲ(2.)*fmax);
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n) * fmax * ਵਰਗਮੂਲ (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  ਜੇ (incx < 0) 
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX    
    x -= (n - 1) * incx * 2;
#ਨਹੀਂਤਾਂ
    x -= (n - 1) * incx;
#ਜੇ_ਅੰਤ
  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE, n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;
}

#ਨਹੀਂਤਾਂ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX
FLOAT CNAME(blasint n, ਖਾਲੀ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ਨਹੀਂਤਾਂ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#ਜੇ_ਅੰਤ

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ਜੇ (n <= 0) ਵਾਪਸ 0.;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (n == 1)
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
    ਵਾਪਸ fabs(x[0]);
#ਨਹੀਂਤਾਂ
    ਵਾਪਸ fabsf(x[0]);
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  ਜੇ (incx == 0)
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n)*fabs(x[0]));
#ਨਹੀਂਤਾਂ
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n)*fabsf(x[0]));
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  {
  ਦੁੱਗਣਾ fr=fabs(x[0]);
  ਦੁੱਗਣਾ fi=fabs(x[1]);
  ਦੁੱਗਣਾ fmin=MIN(fr,fi);
  ਦੁੱਗਣਾ fmax=MAX(fr,fi);
  ਜੇ (fmax==0.) ਵਾਪਸ(fmax);
  ਜੇ (fmax==fmin) ਵਾਪਸ(ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n)*ਵਰਗਮੂਲ(2.)*fmax);
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦੁੱਗਣਾ)n) * fmax * ਵਰਗਮੂਲ (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ਨਹੀਂਤਾਂ
  {
  ਦਸ਼ਮਲਵ fr=fabs(x[0]);
  ਦਸ਼ਮਲਵ fi=fabs(x[1]);
  ਦਸ਼ਮਲਵ fmin=MIN(fr,fi);
  ਦਸ਼ਮਲਵ fmax=MAX(fr,fi);
  ਜੇ (fmax==0.) ਵਾਪਸ(fmax);
  ਜੇ (fmax==fmin) ਵਾਪਸ(ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n)*ਵਰਗਮੂਲ(2.)*fmax);
  ਵਾਪਸ (ਵਰਗਮੂਲ((ਦਸ਼ਮਲਵ)n) * fmax * ਵਰਗਮੂਲ (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  ਜੇ (incx < 0) 
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX    
    x -= (n - 1) * incx * 2;
#ਨਹੀਂਤਾਂ
    x -= (n - 1) * incx;
#ਜੇ_ਅੰਤ
  
  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE, n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;
}

#ਜੇ_ਅੰਤ
