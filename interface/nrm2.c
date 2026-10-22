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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  如果 (n <= 0) 返回 0.;

#ifndef COMPLEX
  如果 (n == 1)
#ifdef DOUBLE
    返回 fabs(x[0]);
#否则
    返回 fabsf(x[0]);
#endif
#endif

  如果 (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  返回 (sqrt((double)n)*fabs(x[0]));
#否则
  返回 (sqrt((float)n)*fabsf(x[0]));
#endif
#否则
#ifdef DOUBLE
  {
  double fr=fabs(x[0]);
  double fi=fabs(x[1]);
  double fmin=MIN(fr,fi);
  double fmax=MAX(fr,fi);
  如果 (fmax==0.) 返回(fmax);
  如果 (fmax==fmin) 返回(sqrt((double)n)*sqrt(2.)*fmax);
  返回 (sqrt((double)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#否则
  {
  float fr=fabs(x[0]);
  float fi=fabs(x[1]);
  float fmin=MIN(fr,fi);
  float fmax=MAX(fr,fi);
  如果 (fmax==0.) 返回(fmax);
  如果 (fmax==fmin) 返回(sqrt((float)n)*sqrt(2.)*fmax);
  返回 (sqrt((float)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  如果 (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#否则
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  返回 ret;
}

#否则

#ifdef COMPLEX
FLOAT CNAME(blasint n, void *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#否则
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  如果 (n <= 0) 返回 0.;

#ifndef COMPLEX
  如果 (n == 1)
#ifdef DOUBLE
    返回 fabs(x[0]);
#否则
    返回 fabsf(x[0]);
#endif
#endif

  如果 (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  返回 (sqrt((double)n)*fabs(x[0]));
#否则
  返回 (sqrt((float)n)*fabsf(x[0]));
#endif
#否则
#ifdef DOUBLE
  {
  double fr=fabs(x[0]);
  double fi=fabs(x[1]);
  double fmin=MIN(fr,fi);
  double fmax=MAX(fr,fi);
  如果 (fmax==0.) 返回(fmax);
  如果 (fmax==fmin) 返回(sqrt((double)n)*sqrt(2.)*fmax);
  返回 (sqrt((double)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#否则
  {
  float fr=fabs(x[0]);
  float fi=fabs(x[1]);
  float fmin=MIN(fr,fi);
  float fmax=MAX(fr,fi);
  如果 (fmax==0.) 返回(fmax);
  如果 (fmax==fmin) 返回(sqrt((float)n)*sqrt(2.)*fmax);
  返回 (sqrt((float)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  如果 (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#否则
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  返回 ret;
}

#endif
