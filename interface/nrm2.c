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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  എങ്കിൽ (n <= 0) തിരികെ 0.;

#ifndef COMPLEX
  എങ്കിൽ (n == 1)
#ifdef DOUBLE
    തിരികെ fabs(x[0]);
#അല്ലെങ്കിൽ
    തിരികെ fabsf(x[0]);
#endif
#endif

  എങ്കിൽ (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  തിരികെ (sqrt((ഇരട്ട)n)*fabs(x[0]));
#അല്ലെങ്കിൽ
  തിരികെ (sqrt((ദശാംശം)n)*fabsf(x[0]));
#endif
#അല്ലെങ്കിൽ
#ifdef DOUBLE
  {
  ഇരട്ട fr=fabs(x[0]);
  ഇരട്ട fi=fabs(x[1]);
  ഇരട്ട fmin=MIN(fr,fi);
  ഇരട്ട fmax=MAX(fr,fi);
  എങ്കിൽ (fmax==0.) തിരികെ(fmax);
  എങ്കിൽ (fmax==fmin) തിരികെ(sqrt((ഇരട്ട)n)*sqrt(2.)*fmax);
  തിരികെ (sqrt((ഇരട്ട)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#അല്ലെങ്കിൽ
  {
  ദശാംശം fr=fabs(x[0]);
  ദശാംശം fi=fabs(x[1]);
  ദശാംശം fmin=MIN(fr,fi);
  ദശാംശം fmax=MAX(fr,fi);
  എങ്കിൽ (fmax==0.) തിരികെ(fmax);
  എങ്കിൽ (fmax==fmin) തിരികെ(sqrt((ദശാംശം)n)*sqrt(2.)*fmax);
  തിരികെ (sqrt((ദശാംശം)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  എങ്കിൽ (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#അല്ലെങ്കിൽ
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  തിരികെ ret;
}

#അല്ലെങ്കിൽ

#ifdef COMPLEX
FLOAT CNAME(blasint n, ശൂന്യം *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#അല്ലെങ്കിൽ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  എങ്കിൽ (n <= 0) തിരികെ 0.;

#ifndef COMPLEX
  എങ്കിൽ (n == 1)
#ifdef DOUBLE
    തിരികെ fabs(x[0]);
#അല്ലെങ്കിൽ
    തിരികെ fabsf(x[0]);
#endif
#endif

  എങ്കിൽ (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  തിരികെ (sqrt((ഇരട്ട)n)*fabs(x[0]));
#അല്ലെങ്കിൽ
  തിരികെ (sqrt((ദശാംശം)n)*fabsf(x[0]));
#endif
#അല്ലെങ്കിൽ
#ifdef DOUBLE
  {
  ഇരട്ട fr=fabs(x[0]);
  ഇരട്ട fi=fabs(x[1]);
  ഇരട്ട fmin=MIN(fr,fi);
  ഇരട്ട fmax=MAX(fr,fi);
  എങ്കിൽ (fmax==0.) തിരികെ(fmax);
  എങ്കിൽ (fmax==fmin) തിരികെ(sqrt((ഇരട്ട)n)*sqrt(2.)*fmax);
  തിരികെ (sqrt((ഇരട്ട)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#അല്ലെങ്കിൽ
  {
  ദശാംശം fr=fabs(x[0]);
  ദശാംശം fi=fabs(x[1]);
  ദശാംശം fmin=MIN(fr,fi);
  ദശാംശം fmax=MAX(fr,fi);
  എങ്കിൽ (fmax==0.) തിരികെ(fmax);
  എങ്കിൽ (fmax==fmin) തിരികെ(sqrt((ദശാംശം)n)*sqrt(2.)*fmax);
  തിരികെ (sqrt((ദശാംശം)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  എങ്കിൽ (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#അല്ലെങ്കിൽ
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  തിരികെ ret;
}

#endif
