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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XSYMV "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZSYMV "
#وإلا
#تعریف_کړه ERROR_NAME "CSYMV "
#د_که_پای

OPENBLAS_EXPORT

فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  صحيح (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#که_تعریف_شوی SMP
  صحيح (*symv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#د_که_پای

  blasint info;
  صحيح uplo;
  FLOAT *بفر;
#که_تعریف_شوی SMP
  صحيح nthreads;
#د_که_پای

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (incy == 0)          info = 10;
  إذا (incx == 0)          info =  7;
  إذا (lda  < MAX(1, n))   info =  5;
  إذا (n < 0)              info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (n == 0) أرجع;

  إذا ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, مطلق_قیمت(incy), عدم, 0, عدم, 0);

  إذا ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  إذا (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

  بفر = (FLOAT *)blas_حافظه_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#که_تعریف_شوی SMP
  nthreads = عدد_cpu_avail(2);

  إذا (n < 50) nthreads = 1;
  إذا (nthreads > 2 && n < 500) nthreads = 2;

  إذا (nthreads == 1) {
#د_که_پای

  (symv[uplo])(n, n, alpha_r, alpha_i, a, lda, b, incx, c, incy, بفر);

#که_تعریف_شوی SMP
  } وإلا {

    (symv_thread[uplo])(n, ALPHA, a, lda, b, incx, c, incy, بفر, nthreads);

  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(4, n * n / 2 + 2 * n,  2 * n * n);

  IDEBUG_پای_ټکی;

  أرجع;
}
