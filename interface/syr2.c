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
#నిర్వచించు ERROR_NAME "QSYR2 "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DSYR2 "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SSYR2 "
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*syr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#నిర్వచితమైతే XDOUBLE
  qsyr2_U, qsyr2_L,
#లేదావేళ defined(DOUBLE)
  dsyr2_U, dsyr2_L,
#లేకపోతే
  ssyr2_U, ssyr2_L,
#అయితే_ముగింపు
};

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*syr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  qsyr2_thread_U, qsyr2_thread_L,
#లేదావేళ defined(DOUBLE)
  dsyr2_thread_U, dsyr2_thread_L,
#లేకపోతే
  ssyr2_thread_U, ssyr2_thread_L,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a, blasint *LDA){

  అక్షరం uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

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

  ఒకవేళ (lda  < MAX(1, n))   info =  9;
  ఒకవేళ (incy == 0)          info =  7;
  ఒకవేళ (incx == 0)          info =  5;
  ఒకవేళ (n < 0)              info =  2;
  ఒకవేళ (uplo  < 0)          info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy, FLOAT *a, blasint lda) {

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

    ఒకవేళ (lda  < MAX(1, n))   info =  9;
    ఒకవేళ (incy == 0)          info =  7;
    ఒకవేళ (incx == 0)          info =  5;
    ఒకవేళ (n < 0)              info =  2;
    ఒకవేళ (uplo  < 0)          info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {

    ఒకవేళ (Uplo == Cblas_పై) uplo  = 1;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 0;

    info = -1;

    ఒకవేళ (lda  < MAX(1, n))   info =  9;
    ఒకవేళ (incy == 0)          info =  7;
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

  ఒకవేళ (incx == 1 && incy == 1 && n < 100) {
    blasint i;
    ఒకవేళ (!uplo) {
      కొరకు (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, ఖాళీ, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, ఖాళీ, 0);
        a += lda;
      }
    } లేకపోతే {
      కొరకు (i = 0; i < n; i++){
        AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, ఖాళీ, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, ఖాళీ, 0);
        a += 1 + lda;
      }
    }
    తిరిగి;
  }

	  
  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0 ) x -= (n - 1) * incx;
  ఒకవేళ (incy < 0 ) y -= (n - 1) * incy;

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#నిర్వచితమైతే SMP
  nthreads = అంకె_cpu_avail(2);

  ఒకవేళ (n < 50) nthreads = 1;
  ఒకవేళ (nthreads > 2 && n < 500) nthreads = 2;

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

    (syr2[uplo])(n, alpha, x, incx, y, incy, a, lda, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (syr2_thread[uplo])(n, alpha, x, incx, y, incy, a, lda, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_చివర;

  తిరిగి;
}
