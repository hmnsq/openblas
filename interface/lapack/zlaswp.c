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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

ساكن صحيح (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#اگر_تعریف_شده XDOUBLE
  xlaswp_plus, xlaswp_minus,
#وإلاإذا defined(DOUBLE)
  zlaswp_plus, zlaswp_minus,
#وإلا
  claswp_plus, claswp_minus,
#پایان_اگر
};

OPENBLAS_EXPORT
فراغ NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  صحيح flag;

#اگر_تعریف_شده SMP
  صحيح نما;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  صحيح nthreads;
#پایان_اگر

  PRINT_DEBUG_NAME;

  إذا (incx == 0 || n <= 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  flag = (incx < 0);

#اگر_تعریف_شده SMP
  إذا (n < 64) {
    nthreads = 1;
  } وإلا {
    nthreads = رقم_cpu_avail(2);
  }

  إذا (nthreads == 1) {
#پایان_اگر

  (laswp[flag])(n, k1, k2, ZERO, ZERO, a, lda, عدم, 0, ipiv, incx);

#اگر_تعریف_شده SMP
  } وإلا {

#اگر_تعریف_شده XDOUBLE
  نما  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  نما  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  نما  =  BLAS_SINGLE  | BLAS_COMPLEX;
#پایان_اگر

  blas_level1_thread(نما, n, k1, k2, dummyalpha, a, lda, عدم, 0, ipiv, incx, (صحيح(*)(فراغ))laswp[flag], nthreads);
  }
#پایان_اگر

  FUNCTION_PROFILE_انتها(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_انتها;

  أرجع;
}
