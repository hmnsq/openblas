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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  যদি (n <= 0) ঘূৰাই 0.;

#ifndef COMPLEX
  যদি (n == 1)
#ifdef DOUBLE
    ঘূৰাই fabs(x[0]);
#অন্যথা
    ঘূৰাই fabsf(x[0]);
#endif
#endif

  যদি (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ঘূৰাই (sqrt((দ্ৱিগুণ)n)*fabs(x[0]));
#অন্যথা
  ঘূৰাই (sqrt((ভিন্ন)n)*fabsf(x[0]));
#endif
#অন্যথা
#ifdef DOUBLE
  {
  দ্ৱিগুণ fr=fabs(x[0]);
  দ্ৱিগুণ fi=fabs(x[1]);
  দ্ৱিগুণ fmin=MIN(fr,fi);
  দ্ৱিগুণ fmax=MAX(fr,fi);
  যদি (fmax==0.) ঘূৰাই(fmax);
  যদি (fmax==fmin) ঘূৰাই(sqrt((দ্ৱিগুণ)n)*sqrt(2.)*fmax);
  ঘূৰাই (sqrt((দ্ৱিগুণ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#অন্যথা
  {
  ভিন্ন fr=fabs(x[0]);
  ভিন্ন fi=fabs(x[1]);
  ভিন্ন fmin=MIN(fr,fi);
  ভিন্ন fmax=MAX(fr,fi);
  যদি (fmax==0.) ঘূৰাই(fmax);
  যদি (fmax==fmin) ঘূৰাই(sqrt((ভিন্ন)n)*sqrt(2.)*fmax);
  ঘূৰাই (sqrt((ভিন্ন)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  যদি (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#অন্যথা
    x -= (n - 1) * incx;
#endif
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ঘূৰাই ret;
}

#অন্যথা

#ifdef COMPLEX
FLOAT CNAME(blasint n, শূন্য *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#অন্যথা
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  যদি (n <= 0) ঘূৰাই 0.;

#ifndef COMPLEX
  যদি (n == 1)
#ifdef DOUBLE
    ঘূৰাই fabs(x[0]);
#অন্যথা
    ঘূৰাই fabsf(x[0]);
#endif
#endif

  যদি (incx == 0)
#ifndef COMPLEX
#ifdef DOUBLE
  ঘূৰাই (sqrt((দ্ৱিগুণ)n)*fabs(x[0]));
#অন্যথা
  ঘূৰাই (sqrt((ভিন্ন)n)*fabsf(x[0]));
#endif
#অন্যথা
#ifdef DOUBLE
  {
  দ্ৱিগুণ fr=fabs(x[0]);
  দ্ৱিগুণ fi=fabs(x[1]);
  দ্ৱিগুণ fmin=MIN(fr,fi);
  দ্ৱিগুণ fmax=MAX(fr,fi);
  যদি (fmax==0.) ঘূৰাই(fmax);
  যদি (fmax==fmin) ঘূৰাই(sqrt((দ্ৱিগুণ)n)*sqrt(2.)*fmax);
  ঘূৰাই (sqrt((দ্ৱিগুণ)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#অন্যথা
  {
  ভিন্ন fr=fabs(x[0]);
  ভিন্ন fi=fabs(x[1]);
  ভিন্ন fmin=MIN(fr,fi);
  ভিন্ন fmax=MAX(fr,fi);
  যদি (fmax==0.) ঘূৰাই(fmax);
  যদি (fmax==fmin) ঘূৰাই(sqrt((ভিন্ন)n)*sqrt(2.)*fmax);
  ঘূৰাই (sqrt((ভিন্ন)n) * fmax * sqrt (1. + (fmin/fmax)*(fmin/fmax)));
  }
#endif
#endif

  যদি (incx < 0) 
#ifdef COMPLEX    
    x -= (n - 1) * incx * 2;
#অন্যথা
    x -= (n - 1) * incx;
#endif
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 2 * n);

  IDEBUG_END;

  ঘূৰাই ret;
}

#endif
