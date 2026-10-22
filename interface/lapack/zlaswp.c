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

ساكن صحيح (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#که_تعریف_شوی XDOUBLE
  xlaswp_plus, xlaswp_minus,
#وإلاإذا defined(DOUBLE)
  zlaswp_plus, zlaswp_minus,
#وإلا
  claswp_plus, claswp_minus,
#د_که_پای
};

OPENBLAS_EXPORT
فراغ NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  صحيح flag;

#که_تعریف_شوی SMP
  صحيح ډېرپېښ;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_NAME;

  إذا (incx == 0 || n <= 0) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  flag = (incx < 0);

#که_تعریف_شوی SMP
  إذا (n < 64) {
    nthreads = 1;
  } وإلا {
    nthreads = عدد_cpu_avail(2);
  }

  إذا (nthreads == 1) {
#د_که_پای

  (laswp[flag])(n, k1, k2, ZERO, ZERO, a, lda, عدم, 0, ipiv, incx);

#که_تعریف_شوی SMP
  } وإلا {

#که_تعریف_شوی XDOUBLE
  ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای

  blas_level1_thread(ډېرپېښ, n, k1, k2, dummyalpha, a, lda, عدم, 0, ipiv, incx, (صحيح(*)(فراغ))laswp[flag], nthreads);
  }
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_پای_ټکی;

  أرجع;
}
