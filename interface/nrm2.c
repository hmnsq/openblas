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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع 0.;

#که_نه_تعریف_شوی COMPLEX
  إذا (n == 1)
#که_تعریف_شوی DOUBLE
    أرجع fabs(x[0]);
#وإلا
    أرجع fabsf(x[0]);
#د_که_پای
#د_که_پای

  إذا (incx == 0)
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی DOUBLE
  أرجع (جذر((مزدوج)n)*fabs(x[0]));
#وإلا
  أرجع (جذر((عائم)n)*fabsf(x[0]));
#د_که_پای
#وإلا
#که_تعریف_شوی DOUBLE
  {
  مزدوج fr=fabs(x[0]);
  مزدوج fi=fabs(x[1]);
  مزدوج fmin=MIN(fr,fi);
  مزدوج fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(جذر((مزدوج)n)*جذر(2.)*fmax);
  أرجع (جذر((مزدوج)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وإلا
  {
  عائم fr=fabs(x[0]);
  عائم fi=fabs(x[1]);
  عائم fmin=MIN(fr,fi);
  عائم fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(جذر((عائم)n)*جذر(2.)*fmax);
  أرجع (جذر((عائم)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#د_که_پای
#د_که_پای

  إذا (incx < 0) 
#که_تعریف_شوی COMPLEX    
    x -= (n - 1) * incx * 2;
#وإلا
    x -= (n - 1) * incx;
#د_که_پای
  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE, n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع ret;
}

#وإلا

#که_تعریف_شوی COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#د_که_پای

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0.;

#که_نه_تعریف_شوی COMPLEX
  إذا (n == 1)
#که_تعریف_شوی DOUBLE
    أرجع fabs(x[0]);
#وإلا
    أرجع fabsf(x[0]);
#د_که_پای
#د_که_پای

  إذا (incx == 0)
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی DOUBLE
  أرجع (جذر((مزدوج)n)*fabs(x[0]));
#وإلا
  أرجع (جذر((عائم)n)*fabsf(x[0]));
#د_که_پای
#وإلا
#که_تعریف_شوی DOUBLE
  {
  مزدوج fr=fabs(x[0]);
  مزدوج fi=fabs(x[1]);
  مزدوج fmin=MIN(fr,fi);
  مزدوج fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(جذر((مزدوج)n)*جذر(2.)*fmax);
  أرجع (جذر((مزدوج)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#وإلا
  {
  عائم fr=fabs(x[0]);
  عائم fi=fabs(x[1]);
  عائم fmin=MIN(fr,fi);
  عائم fmax=MAX(fr,fi);
  إذا (fmax==0.) أرجع(fmax);
  إذا (fmax==fmin) أرجع(جذر((عائم)n)*جذر(2.)*fmax);
  أرجع (جذر((عائم)n) * fmax * جذر (1. + (fmin/fmax)*(fmin/fmax)));
  }
#د_که_پای
#د_که_پای

  إذا (incx < 0) 
#که_تعریف_شوی COMPLEX    
    x -= (n - 1) * incx * 2;
#وإلا
    x -= (n - 1) * incx;
#د_که_پای
  
  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE, n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع ret;
}

#د_که_پای
