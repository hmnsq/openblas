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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  अगर (n <= 0) वापस 0.;

#अगर_अघोषित COMPLEX
  अगर (n == 1)
#अगर_घोषित DOUBLE
    वापस भ_असल(x[0]);
#अन्यथा
    वापस fabsf(x[0]);
#पूर्ण_अगर
#पूर्ण_अगर

  अगर (incx == 0)
#अगर_अघोषित COMPLEX
#अगर_घोषित DOUBLE
  वापस (वर्ग_मूल((द्विगुन)n)*भ_असल(x[0]));
#अन्यथा
  वापस (वर्ग_मूल((भग्न)n)*fabsf(x[0]));
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित DOUBLE
  {
  द्विगुन fr=भ_असल(x[0]);
  द्विगुन fi=भ_असल(x[1]);
  द्विगुन fmin=MIN(fr,fi);
  द्विगुन fmax=MAX(fr,fi);
  अगर (fmax==0.) वापस(fmax);
  अगर (fmax==fmin) वापस(वर्ग_मूल((द्विगुन)n)*वर्ग_मूल(2.)*fmax);
  वापस (वर्ग_मूल((द्विगुन)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#अन्यथा
  {
  भग्न fr=भ_असल(x[0]);
  भग्न fi=भ_असल(x[1]);
  भग्न fmin=MIN(fr,fi);
  भग्न fmax=MAX(fr,fi);
  अगर (fmax==0.) वापस(fmax);
  अगर (fmax==fmin) वापस(वर्ग_मूल((भग्न)n)*वर्ग_मूल(2.)*fmax);
  वापस (वर्ग_मूल((भग्न)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#पूर्ण_अगर
#पूर्ण_अगर

  अगर (incx < 0) 
#अगर_घोषित COMPLEX    
    x -= (n - 1) * incx * 2;
#अन्यथा
    x -= (n - 1) * incx;
#पूर्ण_अगर
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  वापस ret;
}

#अन्यथा

#अगर_घोषित COMPLEX
FLOAT CNAME(blasint n, व्योम *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#पूर्ण_अगर

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  अगर (n <= 0) वापस 0.;

#अगर_अघोषित COMPLEX
  अगर (n == 1)
#अगर_घोषित DOUBLE
    वापस भ_असल(x[0]);
#अन्यथा
    वापस fabsf(x[0]);
#पूर्ण_अगर
#पूर्ण_अगर

  अगर (incx == 0)
#अगर_अघोषित COMPLEX
#अगर_घोषित DOUBLE
  वापस (वर्ग_मूल((द्विगुन)n)*भ_असल(x[0]));
#अन्यथा
  वापस (वर्ग_मूल((भग्न)n)*fabsf(x[0]));
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित DOUBLE
  {
  द्विगुन fr=भ_असल(x[0]);
  द्विगुन fi=भ_असल(x[1]);
  द्विगुन fmin=MIN(fr,fi);
  द्विगुन fmax=MAX(fr,fi);
  अगर (fmax==0.) वापस(fmax);
  अगर (fmax==fmin) वापस(वर्ग_मूल((द्विगुन)n)*वर्ग_मूल(2.)*fmax);
  वापस (वर्ग_मूल((द्विगुन)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#अन्यथा
  {
  भग्न fr=भ_असल(x[0]);
  भग्न fi=भ_असल(x[1]);
  भग्न fmin=MIN(fr,fi);
  भग्न fmax=MAX(fr,fi);
  अगर (fmax==0.) वापस(fmax);
  अगर (fmax==fmin) वापस(वर्ग_मूल((भग्न)n)*वर्ग_मूल(2.)*fmax);
  वापस (वर्ग_मूल((भग्न)n) * fmax * वर्ग_मूल (1. + (fmin/fmax)*(fmin/fmax)));
  }
#पूर्ण_अगर
#पूर्ण_अगर

  अगर (incx < 0) 
#अगर_घोषित COMPLEX    
    x -= (n - 1) * incx * 2;
#अन्यथा
    x -= (n - 1) * incx;
#पूर्ण_अगर
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  वापस ret;
}

#पूर्ण_अगर
