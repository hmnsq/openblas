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
#అన్తర్భావ <ctype.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XHER  "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZHER  "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CHER  "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#నిర్వచితమైతే XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#లేదావేళ defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#లేకపోతే
  cher_U, cher_L, cher_V, cher_M,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#లేదావేళ defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#లేకపోతే
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  అక్షరం uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  పూర్ణాంకం uplo;
  FLOAT *బఫరు;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ఒకవేళ (lda  < MAX(1, n))   info =  7;
  ఒకవేళ (incx == 0)          info =  5;
  ఒకవేళ (n < 0)              info =  2;
  ఒకవేళ (uplo  < 0)          info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo, blasint n, FLOAT alpha, శూన్య *vx, blasint incx, శూన్య *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *బఫరు;
  పూర్ణాంకం uplo;
  blasint info;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {

    ఒకవేళ (Uplo == Cblas_పై) uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 1;

    info = -1;

    ఒకవేళ (lda  < MAX(1, n))   info =  7;
    ఒకవేళ (incx == 0)          info =  5;
    ఒకవేళ (n < 0)              info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;

  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {

    ఒకవేళ (Uplo == Cblas_పై) uplo  = 3;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 2;

    info = -1;

    ఒకవేళ (lda  < MAX(1, n))   info =  7;
    ఒకవేళ (incx == 0)          info =  5;
    ఒకవేళ (n < 0)              info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  ఒకవేళ (n == 0) తిరిగి;

  ఒకవేళ (alpha == ZERO) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0 ) x -= (n - 1) * incx * 2;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#నిర్వచితమైతే SMP
  ఒకవేళ (n < 100)
    nthreads = 1;
  లేకపోతే
    nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  (her[uplo])(n, alpha, x, incx, a, lda, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(4, n * n / 2 + n, n * n);

  IDEBUG_చివర;

  తిరిగి;
}
