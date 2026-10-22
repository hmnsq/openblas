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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  اطبع_DEBUG_NAME;

  إذا (n <= 0) أرجع 0.;

#ifndef COMPLEX
  إذا (n == 1)
#ifdef DOUBLE
    أرجع fabs(x[0]);
#وإلا
    أرجع fabsf(x[0]);
#endif
#endif

  إذا (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  أرجع (sqrt((مزدوج)n)*fabs(x[0]));
#وإلا
  أرجع (sqrt((عائم)n)*fabsf(x[0]));
#endif
#وإلا
#ifdef DOUBLE
  {
  مزدوج fr=fabs(x[0]);
  مزدوج fi=fabs(x[1]);
  مزدوج fmin=MIN(fr,fi);
  مزدوج fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(sqrt((مزدوج)n)*sqrt(2.)*fmax);
  أرجع (sqrt((مزدوج)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وإلا
  {
  عائم fr=fabs(x[0]);
  عائم fi=fabs(x[1]);
  عائم fmin=MIN(fr,fi);
  عائم fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(sqrt((عائم)n)*sqrt(2.)*fmax);
  أرجع (sqrt((عائم)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  إذا (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#وإلا
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  أرجع ret;
}

#وإلا

#ifdef COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  اطبع_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0.;

#ifndef COMPLEX
  إذا (n == 1)
#ifdef DOUBLE
    أرجع fabs(x[0]);
#وإلا
    أرجع fabsf(x[0]);
#endif
#endif

  إذا (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  أرجع (sqrt((مزدوج)n)*fabs(x[0]));
#وإلا
  أرجع (sqrt((عائم)n)*fabsf(x[0]));
#endif
#وإلا
#ifdef DOUBLE
  {
  مزدوج fr=fabs(x[0]);
  مزدوج fi=fabs(x[1]);
  مزدوج fmin=MIN(fr,fi);
  مزدوج fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(sqrt((مزدوج)n)*sqrt(2.)*fmax);
  أرجع (sqrt((مزدوج)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وإلا
  {
  عائم fr=fabs(x[0]);
  عائم fi=fabs(x[1]);
  عائم fmin=MIN(fr,fi);
  عائم fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(sqrt((عائم)n)*sqrt(2.)*fmax);
  أرجع (sqrt((عائم)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  إذا (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#وإلا
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  أرجع ret;
}

#endif
