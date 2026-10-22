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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#ifdef XDOUBLE
  xlaswp_plus, xlaswp_minus,
#ନଚେତଯଦି defined(DOUBLE)
  zlaswp_plus, zlaswp_minus,
#ଅନ୍ୟଥା
  claswp_plus, claswp_minus,
#endif
};

OPENBLAS_EXPORT
ଶୂନ୍ୟ NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା flag;

#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  PRINT_DEBUG_NAME;

  ଯଦି (incx == 0 || n <= 0) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  flag = (incx < 0);

#ifdef SMP
  ଯଦି (n < 64) {
    nthreads = 1;
  } ଅନ୍ୟଥା {
    nthreads = ଅଙ୍କ_cpu_avail(2);
  }

  ଯଦି (nthreads == 1) {
#endif

  (laswp[flag])(n, k1, k2, ZERO, ZERO, a, lda, ଖାଲି, 0, ipiv, incx);

#ifdef SMP
  } ଅନ୍ୟଥା {

#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif

  blas_level1_thread(mode, n, k1, k2, dummyalpha, a, lda, ଖାଲି, 0, ipiv, incx, (ପୂର୍ଣ୍ଣସଂଖ୍ୟା(*)(ଶୂନ୍ୟ))laswp[flag], nthreads);
  }
#endif

  FUNCTION_PROFILE_END(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
