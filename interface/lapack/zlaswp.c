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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#నిర్వచితమైతే XDOUBLE
  xlaswp_plus, xlaswp_minus,
#లేదావేళ defined(DOUBLE)
  zlaswp_plus, zlaswp_minus,
#లేకపోతే
  claswp_plus, claswp_minus,
#అయితే_ముగింపు
};

OPENBLAS_EXPORT
శూన్య NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  పూర్ణాంకం flag;

#నిర్వచితమైతే SMP
  పూర్ణాంకం బహుళకం;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  ఒకవేళ (incx == 0 || n <= 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  flag = (incx < 0);

#నిర్వచితమైతే SMP
  ఒకవేళ (n < 64) {
    nthreads = 1;
  } లేకపోతే {
    nthreads = అంకె_cpu_avail(2);
  }

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (laswp[flag])(n, k1, k2, ZERO, ZERO, a, lda, ఖాళీ, 0, ipiv, incx);

#నిర్వచితమైతే SMP
  } లేకపోతే {

#నిర్వచితమైతే XDOUBLE
  బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు

  blas_level1_thread(బహుళకం, n, k1, k2, dummyalpha, a, lda, ఖాళీ, 0, ipiv, incx, (పూర్ణాంకం(*)(శూన్య))laswp[flag], nthreads);
  }
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_చివర;

  తిరిగి;
}
