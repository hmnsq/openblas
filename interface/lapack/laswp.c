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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qlaswp_plus, qlaswp_minus,
#ਜਾਂਜੇ defined(DOUBLE)
  dlaswp_plus, dlaswp_minus,
#ਨਹੀਂਤਾਂ
  slaswp_plus, slaswp_minus,
#ਜੇ_ਅੰਤ
};

OPENBLAS_EXPORT
ਖਾਲੀ NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  ਪੂਰਨ flag;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ ਬਹੁਲਕ, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  ਜੇ (incx == 0 || n <= 0) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  flag = (incx < 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਜੇ (n < 64) {
    nthreads = 1;
  } ਨਹੀਂਤਾਂ {
    nthreads = ਅੰਕ_cpu_avail(1);
  }

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  (laswp[flag])(n, k1, k2, ZERO, a, lda, ਨਿਰਾ, 0, ipiv, incx);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ

  blas_level1_thread(ਬਹੁਲਕ, n, k1, k2, dummyalpha,
		     a, lda, ਨਿਰਾ, 0, ipiv, incx,
		     (ਪੂਰਨ(*)(ਖਾਲੀ))laswp[flag], nthreads);
  }
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;

}
