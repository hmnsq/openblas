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

अटल पूर्णांक (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#अगर_घोषित XDOUBLE
  qlaswp_plus, qlaswp_minus,
#या_अगर defined(DOUBLE)
  dlaswp_plus, dlaswp_minus,
#अन्यथा
  slaswp_plus, slaswp_minus,
#पूर्ण_अगर
};

OPENBLAS_EXPORT
व्योम NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  पूर्णांक flag;

#अगर_घोषित SMP
  पूर्णांक mode, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  अगर (incx == 0 || n <= 0) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  flag = (incx < 0);

#अगर_घोषित SMP
  अगर (n < 64) {
    nthreads = 1;
  } अन्यथा {
    nthreads = अंक_cpu_avail(1);
  }

  अगर (nthreads == 1) {
#पूर्ण_अगर

  (laswp[flag])(n, k1, k2, ZERO, a, lda, शून्य, 0, ipiv, incx);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_घोषित XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर

  blas_level1_thread(mode, n, k1, k2, dummyalpha,
		     a, lda, शून्य, 0, ipiv, incx,
		     (पूर्णांक(*)(व्योम))laswp[flag], nthreads);
  }
#पूर्ण_अगर

  FUNCTION_PROFILE_END(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_END;

  वापस;

}
