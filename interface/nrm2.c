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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ଯଦି (n <= 0) ଫେରନ୍ତୁ 0.;

#ifndef COMPLEX
  ଯଦି (n == 1)
#ifdef DOUBLE
    ଫେରନ୍ତୁ fabs(x[0]);
#ଅନ୍ୟଥା
    ଫେରନ୍ତୁ fabsf(x[0]);
#endif
#endif

  ଯଦି (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ଫେରନ୍ତୁ (sqrt((ଦ୍ୱିଗୁଣ)n)*fabs(x[0]));
#ଅନ୍ୟଥା
  ଫେରନ୍ତୁ (sqrt((ଦଶମିକ)n)*fabsf(x[0]));
#endif
#ଅନ୍ୟଥା
#ifdef DOUBLE
  {
  ଦ୍ୱିଗୁଣ fr=fabs(x[0]);
  ଦ୍ୱିଗୁଣ fi=fabs(x[1]);
  ଦ୍ୱିଗୁଣ fmin=MIN(fr,fi);
  ଦ୍ୱିଗୁଣ fmax=MAX(fr,fi);
  ଯଦି (fmax==0.) ଫେରନ୍ତୁ(fmax);
  ଯଦି (fmax==fmin) ଫେରନ୍ତୁ(sqrt((ଦ୍ୱିଗୁଣ)n)*sqrt(2.)*fmax);
  ଫେରନ୍ତୁ (sqrt((ଦ୍ୱିଗୁଣ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ଅନ୍ୟଥା
  {
  ଦଶମିକ fr=fabs(x[0]);
  ଦଶମିକ fi=fabs(x[1]);
  ଦଶମିକ fmin=MIN(fr,fi);
  ଦଶମିକ fmax=MAX(fr,fi);
  ଯଦି (fmax==0.) ଫେରନ୍ତୁ(fmax);
  ଯଦି (fmax==fmin) ଫେରନ୍ତୁ(sqrt((ଦଶମିକ)n)*sqrt(2.)*fmax);
  ଫେରନ୍ତୁ (sqrt((ଦଶମିକ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  ଯଦି (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#ଅନ୍ୟଥା
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ ret;
}

#ଅନ୍ୟଥା

#ifdef COMPLEX
FLOAT CNAME(blasint n, ଶୂନ୍ୟ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ଅନ୍ୟଥା
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ଯଦି (n <= 0) ଫେରନ୍ତୁ 0.;

#ifndef COMPLEX
  ଯଦି (n == 1)
#ifdef DOUBLE
    ଫେରନ୍ତୁ fabs(x[0]);
#ଅନ୍ୟଥା
    ଫେରନ୍ତୁ fabsf(x[0]);
#endif
#endif

  ଯଦି (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ଫେରନ୍ତୁ (sqrt((ଦ୍ୱିଗୁଣ)n)*fabs(x[0]));
#ଅନ୍ୟଥା
  ଫେରନ୍ତୁ (sqrt((ଦଶମିକ)n)*fabsf(x[0]));
#endif
#ଅନ୍ୟଥା
#ifdef DOUBLE
  {
  ଦ୍ୱିଗୁଣ fr=fabs(x[0]);
  ଦ୍ୱିଗୁଣ fi=fabs(x[1]);
  ଦ୍ୱିଗୁଣ fmin=MIN(fr,fi);
  ଦ୍ୱିଗୁଣ fmax=MAX(fr,fi);
  ଯଦି (fmax==0.) ଫେରନ୍ତୁ(fmax);
  ଯଦି (fmax==fmin) ଫେରନ୍ତୁ(sqrt((ଦ୍ୱିଗୁଣ)n)*sqrt(2.)*fmax);
  ଫେରନ୍ତୁ (sqrt((ଦ୍ୱିଗୁଣ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ଅନ୍ୟଥା
  {
  ଦଶମିକ fr=fabs(x[0]);
  ଦଶମିକ fi=fabs(x[1]);
  ଦଶମିକ fmin=MIN(fr,fi);
  ଦଶମିକ fmax=MAX(fr,fi);
  ଯଦି (fmax==0.) ଫେରନ୍ତୁ(fmax);
  ଯଦି (fmax==fmin) ଫେରନ୍ତୁ(sqrt((ଦଶମିକ)n)*sqrt(2.)*fmax);
  ଫେରନ୍ତୁ (sqrt((ଦଶମିକ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  ଯଦି (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#ଅନ୍ୟଥା
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ ret;
}

#endif
