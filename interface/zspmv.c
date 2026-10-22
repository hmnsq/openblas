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
#define ERROR_NAME "ZSPMV "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZSPMV "
#অন্যথা
#define ERROR_NAME "CSPMV "
#endif

স্থির  পূর্ণসংখ্যা (*spmv[])(BLASLONG, FLOAT, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
  xspmv_U, xspmv_L,
#নতুবা_যদি defined(DOUBLE)
  zspmv_U, zspmv_L,
#অন্যথা
  cspmv_U, cspmv_L,
#endif
};

#ifdef SMPTEST
স্থির  পূর্ণসংখ্যা (*spmv_thread[])(BLASLONG, FLOAT *, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণসংখ্যা) = {
#ifdef XDOUBLE
  xspmv_thread_U, xspmv_thread_L,
#নতুবা_যদি defined(DOUBLE)
  zspmv_thread_U, zspmv_thread_L,
#অন্যথা
  cspmv_thread_U, cspmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

শূন্য NAME(অক্ষর *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  অক্ষর uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  পূর্ণসংখ্যা uplo;
  FLOAT *বাফার;
#ifdef SMPTEST
  পূর্ণসংখ্যা nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  যদি (incy == 0)          info =  9;
  যদি (incx == 0)          info =  6;
  যদি (n < 0)              info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

  যদি (n == 0) ফেরত;

  যদি ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, abs(incy), ফাঁকা, 0, ফাঁকা, 0);

  যদি ((alpha_r == ZERO) && (alpha_i == ZERO)) ফেরত;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  যদি (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

  বাফার = (FLOAT *)blas_memory_alloc(1);
  যদি (!বাফার) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#ifdef SMPTEST
  nthreads = অঙ্ক_cpu_avail(2);

  যদি (nthreads == 1) {
#endif

  (spmv[uplo])(n, alpha_r, alpha_i, a, b, incx, c, incy, বাফার);

#ifdef SMPTEST

  } অন্যথা {

    (spmv_thread[uplo])(n, ALPHA, a, b, incx, c, incy, বাফার, nthreads);

  }
#endif

  blas_memory_free(বাফার);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ফেরত;
}
