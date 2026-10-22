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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XSYR  "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZSYR  "
#অন্যথা
#define ERROR_NAME "CSYR  "
#endif

স্থির পূর্ণসংখ্যা (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr_U, xsyr_L,
#নতুবা_যদি defined(DOUBLE)
  zsyr_U, zsyr_L,
#অন্যথা
  csyr_U, csyr_L,
#endif
};

#ifdef SMP
স্থির পূর্ণসংখ্যা (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণসংখ্যা) = {
#ifdef XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#নতুবা_যদি defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#অন্যথা
  csyr_thread_U, csyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(অক্ষর *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  অক্ষর uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  পূর্ণসংখ্যা uplo;
  FLOAT *বাফার;
#ifdef SMP
  পূর্ণসংখ্যা nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  যদি (lda  < MAX(1, n))   info =  7;
  যদি (incx == 0)          info =  5;
  যদি (n < 0)              info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }


#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo, পূর্ণসংখ্যা n, শূন্য* valpha, FLOAT *x, পূর্ণসংখ্যা incx, FLOAT *a, পূর্ণসংখ্যা lda) {

  FLOAT *বাফার;
  পূর্ণসংখ্যা uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
#ifdef SMP
  পূর্ণসংখ্যা nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  যদি (order == CblasColMajor) {

    যদি (Uplo == CblasUpper) uplo  = 0;
    যদি (Uplo == CblasLower) uplo  = 1;

    info = -1;

    যদি (lda  < MAX(1, n))   info =  7;
    যদি (incx == 0)          info =  5;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;

  }

  যদি (order == CblasRowMajor) {

    যদি (Uplo == CblasUpper) uplo  = 1;
    যদি (Uplo == CblasLower) uplo  = 0;

    info = -1;

    যদি (lda  < MAX(1, n))   info =  7;
    যদি (incx == 0)          info =  5;
    যদি (n < 0)              info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#endif

  যদি (n == 0) ফেরত;

  যদি ((alpha_r == ZERO) && (alpha_i == ZERO)) ফেরত;

  যদি (incx == 1 && n < 50) {
    blasint i;
    যদি (!uplo) {
      জন্য (i = 0; i < n; i++){
        যদি ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, ফাঁকা, 0);
        }
        a += lda * 2;
      }
    } অন্যথা {
      জন্য (i = 0; i < n; i++){
        যদি ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, ফাঁকা, 0);
        }
        a += 2 + lda * 2;
      }
    }
    ফেরত;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx * 2;

  বাফার = (FLOAT *)blas_memory_alloc(1);
  যদি (!বাফার) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#ifdef SMP
  nthreads = অঙ্ক_cpu_avail(2);

  যদি (n < 50) nthreads = 1;
  যদি (nthreads > 2 && n < 500) nthreads = 2;

  যদি (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, বাফার);

#ifdef SMP
  } অন্যথা {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, বাফার, nthreads);

  }
#endif

  blas_memory_free(বাফার);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ফেরত;
}
