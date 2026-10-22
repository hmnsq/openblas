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

স্থির পূর্ণাঙ্ক (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#ifdef XDOUBLE
  qlaswp_plus, qlaswp_minus,
#নোচেদ্যদি defined(DOUBLE)
  dlaswp_plus, dlaswp_minus,
#অন্যথা
  slaswp_plus, slaswp_minus,
#endif
};

OPENBLAS_EXPORT
শূন্য NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  পূর্ণাঙ্ক flag;

#ifdef SMP
  পূর্ণাঙ্ক mode, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#endif

  PRINT_DEBUG_NAME;

  যদি (incx == 0 || n <= 0) ঘূৰাই;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  flag = (incx < 0);

#ifdef SMP
  যদি (n < 64) {
    nthreads = 1;
  } অন্যথা {
    nthreads = অঙ্ক_cpu_avail(1);
  }

  যদি (nthreads == 1) {
#endif

  (laswp[flag])(n, k1, k2, ZERO, a, lda, রিক্ত, 0, ipiv, incx);

#ifdef SMP
  } অন্যথা {

#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif

  blas_level1_thread(mode, n, k1, k2, dummyalpha,
		     a, lda, রিক্ত, 0, ipiv, incx,
		     (পূর্ণাঙ্ক(*)(শূন্য))laswp[flag], nthreads);
  }
#endif

  FUNCTION_PROFILE_END(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_END;

  ঘূৰাই;

}
