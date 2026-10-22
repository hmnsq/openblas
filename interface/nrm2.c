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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0.;

#ifndef COMPLEX
  ಯದಿ (n == 1)
#ifdef DOUBLE
    ಪ್ರತ್ಯಾಗಮ fabs(x[0]);
#ಅನ್ಯಥಾ
    ಪ್ರತ್ಯಾಗಮ fabsf(x[0]);
#endif
#endif

  ಯದಿ (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ಪ್ರತ್ಯಾಗಮ (sqrt((ದ್ವಿಗುಣ)n)*fabs(x[0]));
#ಅನ್ಯಥಾ
  ಪ್ರತ್ಯಾಗಮ (sqrt((ಭಿನ್ನ)n)*fabsf(x[0]));
#endif
#ಅನ್ಯಥಾ
#ifdef DOUBLE
  {
  ದ್ವಿಗುಣ fr=fabs(x[0]);
  ದ್ವಿಗುಣ fi=fabs(x[1]);
  ದ್ವಿಗುಣ fmin=MIN(fr,fi);
  ದ್ವಿಗುಣ fmax=MAX(fr,fi);
  ಯದಿ (fmax==0.) ಪ್ರತ್ಯಾಗಮ(fmax);
  ಯದಿ (fmax==fmin) ಪ್ರತ್ಯಾಗಮ(sqrt((ದ್ವಿಗುಣ)n)*sqrt(2.)*fmax);
  ಪ್ರತ್ಯಾಗಮ (sqrt((ದ್ವಿಗುಣ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ಅನ್ಯಥಾ
  {
  ಭಿನ್ನ fr=fabs(x[0]);
  ಭಿನ್ನ fi=fabs(x[1]);
  ಭಿನ್ನ fmin=MIN(fr,fi);
  ಭಿನ್ನ fmax=MAX(fr,fi);
  ಯದಿ (fmax==0.) ಪ್ರತ್ಯಾಗಮ(fmax);
  ಯದಿ (fmax==fmin) ಪ್ರತ್ಯಾಗಮ(sqrt((ಭಿನ್ನ)n)*sqrt(2.)*fmax);
  ಪ್ರತ್ಯಾಗಮ (sqrt((ಭಿನ್ನ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  ಯದಿ (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#ಅನ್ಯಥಾ
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;
}

#ಅನ್ಯಥಾ

#ifdef COMPLEX
FLOAT CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ಅನ್ಯಥಾ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0.;

#ifndef COMPLEX
  ಯದಿ (n == 1)
#ifdef DOUBLE
    ಪ್ರತ್ಯಾಗಮ fabs(x[0]);
#ಅನ್ಯಥಾ
    ಪ್ರತ್ಯಾಗಮ fabsf(x[0]);
#endif
#endif

  ಯದಿ (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ಪ್ರತ್ಯಾಗಮ (sqrt((ದ್ವಿಗುಣ)n)*fabs(x[0]));
#ಅನ್ಯಥಾ
  ಪ್ರತ್ಯಾಗಮ (sqrt((ಭಿನ್ನ)n)*fabsf(x[0]));
#endif
#ಅನ್ಯಥಾ
#ifdef DOUBLE
  {
  ದ್ವಿಗುಣ fr=fabs(x[0]);
  ದ್ವಿಗುಣ fi=fabs(x[1]);
  ದ್ವಿಗುಣ fmin=MIN(fr,fi);
  ದ್ವಿಗುಣ fmax=MAX(fr,fi);
  ಯದಿ (fmax==0.) ಪ್ರತ್ಯಾಗಮ(fmax);
  ಯದಿ (fmax==fmin) ಪ್ರತ್ಯಾಗಮ(sqrt((ದ್ವಿಗುಣ)n)*sqrt(2.)*fmax);
  ಪ್ರತ್ಯಾಗಮ (sqrt((ದ್ವಿಗುಣ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#ಅನ್ಯಥಾ
  {
  ಭಿನ್ನ fr=fabs(x[0]);
  ಭಿನ್ನ fi=fabs(x[1]);
  ಭಿನ್ನ fmin=MIN(fr,fi);
  ಭಿನ್ನ fmax=MAX(fr,fi);
  ಯದಿ (fmax==0.) ಪ್ರತ್ಯಾಗಮ(fmax);
  ಯದಿ (fmax==fmin) ಪ್ರತ್ಯಾಗಮ(sqrt((ಭಿನ್ನ)n)*sqrt(2.)*fmax);
  ಪ್ರತ್ಯಾಗಮ (sqrt((ಭಿನ್ನ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  ಯದಿ (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#ಅನ್ಯಥಾ
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;
}

#endif
