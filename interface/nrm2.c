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

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  यदि (n <= 0) प्रत्यागम 0.;

#ifndef COMPLEX
  यदि (n == 1)
#ifdef DOUBLE
    प्रत्यागम भ_असल(x[0]);
#अन्यथा
    प्रत्यागम fabsf(x[0]);
#endif
#endif

  यदि (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  प्रत्यागम (वर्ग_मूल((द्विगुण)n)*भ_असल(x[0]));
#अन्यथा
  प्रत्यागम (वर्ग_मूल((भिन्न)n)*fabsf(x[0]));
#endif
#अन्यथा
#ifdef DOUBLE
  {
  द्विगुण fr=भ_असल(x[0]);
  द्विगुण fi=भ_असल(x[1]);
  द्विगुण fmin=MIN(fr,fi);
  द्विगुण fmax=MAX(fr,fi);
  यदि (fmax==0.) प्रत्यागम(fmax);
  यदि (fmax==fmin) प्रत्यागम(वर्ग_मूल((द्विगुण)n)*वर्ग_मूल(2.)*fmax);
  प्रत्यागम (वर्ग_मूल((द्विगुण)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#अन्यथा
  {
  भिन्न fr=भ_असल(x[0]);
  भिन्न fi=भ_असल(x[1]);
  भिन्न fmin=MIN(fr,fi);
  भिन्न fmax=MAX(fr,fi);
  यदि (fmax==0.) प्रत्यागम(fmax);
  यदि (fmax==fmin) प्रत्यागम(वर्ग_मूल((भिन्न)n)*वर्ग_मूल(2.)*fmax);
  प्रत्यागम (वर्ग_मूल((भिन्न)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  यदि (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#अन्यथा
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  प्रत्यागम ret;
}

#अन्यथा

#ifdef COMPLEX
FLOAT CNAME(blasint n, शून्य *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  यदि (n <= 0) प्रत्यागम 0.;

#ifndef COMPLEX
  यदि (n == 1)
#ifdef DOUBLE
    प्रत्यागम भ_असल(x[0]);
#अन्यथा
    प्रत्यागम fabsf(x[0]);
#endif
#endif

  यदि (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  प्रत्यागम (वर्ग_मूल((द्विगुण)n)*भ_असल(x[0]));
#अन्यथा
  प्रत्यागम (वर्ग_मूल((भिन्न)n)*fabsf(x[0]));
#endif
#अन्यथा
#ifdef DOUBLE
  {
  द्विगुण fr=भ_असल(x[0]);
  द्विगुण fi=भ_असल(x[1]);
  द्विगुण fmin=MIN(fr,fi);
  द्विगुण fmax=MAX(fr,fi);
  यदि (fmax==0.) प्रत्यागम(fmax);
  यदि (fmax==fmin) प्रत्यागम(वर्ग_मूल((द्विगुण)n)*वर्ग_मूल(2.)*fmax);
  प्रत्यागम (वर्ग_मूल((द्विगुण)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#अन्यथा
  {
  भिन्न fr=भ_असल(x[0]);
  भिन्न fi=भ_असल(x[1]);
  भिन्न fmin=MIN(fr,fi);
  भिन्न fmax=MAX(fr,fi);
  यदि (fmax==0.) प्रत्यागम(fmax);
  यदि (fmax==fmin) प्रत्यागम(वर्ग_मूल((भिन्न)n)*वर्ग_मूल(2.)*fmax);
  प्रत्यागम (वर्ग_मूल((भिन्न)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  यदि (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#अन्यथा
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  प्रत्यागम ret;
}

#endif
