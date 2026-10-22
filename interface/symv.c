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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <ctype.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QSYMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DSYMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SSYMV "
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  ਅੱਖਰ uplo_arg = *UPLO;
  blasint n	= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  ਪੂਰਨ (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#ਜੇ_ਅੰਤ

  blasint info;
  ਪੂਰਨ uplo;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ਜੇ (uplo_arg  == 'U') uplo  = 0;
  ਜੇ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ਜੇ (incy == 0)          info = 10;
  ਜੇ (incx == 0)          info =  7;
  ਜੇ (lda  < MAX(1, n))   info =  5;
  ਜੇ (n < 0)              info =  2;
  ਜੇ (uplo  < 0)          info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_UPLO Uplo, blasint n, FLOAT alpha,
	   FLOAT *a, blasint lda, FLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy) {

  FLOAT *ਬਫਰ;
  ਪੂਰਨ uplo;
  blasint info;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  ਪੂਰਨ (*symv[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    SYMV_U, SYMV_L,
  };

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ (*symv_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
    SYMV_THREAD_U, SYMV_THREAD_L,
  };
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 1;

    info = -1;

    ਜੇ (incy == 0)          info = 10;
    ਜੇ (incx == 0)          info =  7;
    ਜੇ (lda  < MAX(1, n))   info =  5;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 1;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 0;

    info = -1;

    ਜੇ (incy == 0)          info = 10;
    ਜੇ (incx == 0)          info =  7;
    ਜੇ (lda  < MAX(1, n))   info =  5;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  ਜੇ (n == 0) ਵਾਪਸ;

  ਜੇ (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);

  ਜੇ (alpha == ZERO) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0 ) x -= (n - 1) * incx;
  ਜੇ (incy < 0 ) y -= (n - 1) * incy;

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਜੇ (n <200)
	  nthreads=1;
  ਨਹੀਂਤਾਂ
  nthreads = ਅੰਕ_cpu_avail(2);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  (symv[uplo])(n, n, alpha, a, lda, x, incx, y, incy, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    (symv_thread[uplo])(n, alpha, a, lda, x, incx, y, incy, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(1, n * n / 2 + 2 * n,  2 * n * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
