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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

स्थिर पूर्णाङ्क (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#यदि_परिभाषितम् XDOUBLE
  qlaswp_plus, qlaswp_minus,
#नोचेद्यदि defined(DOUBLE)
  dlaswp_plus, dlaswp_minus,
#अन्यथा
  slaswp_plus, slaswp_minus,
#यदि_अन्तः
};

OPENBLAS_EXPORT
शून्य NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  पूर्णाङ्क flag;

#यदि_परिभाषितम् SMP
  पूर्णाङ्क बहुलकम्, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  यदि (incx == 0 || n <= 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  flag = (incx < 0);

#यदि_परिभाषितम् SMP
  यदि (n < 64) {
    nthreads = 1;
  } अन्यथा {
    nthreads = अङ्कः_cpu_avail(1);
  }

  यदि (nthreads == 1) {
#यदि_अन्तः

  (laswp[flag])(n, k1, k2, ZERO, a, lda, रिक्त, 0, ipiv, incx);

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_परिभाषितम् XDOUBLE
  बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः

  blas_level1_thread(बहुलकम्, n, k1, k2, dummyalpha,
		     a, lda, रिक्त, 0, ipiv, incx,
		     (पूर्णाङ्क(*)(शून्य))laswp[flag], nthreads);
  }
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
