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
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QSYR  "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DSYR  "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SSYR  "
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*syr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qsyr_U, qsyr_L,
#ਜਾਂਜੇ defined(DOUBLE)
  dsyr_U, dsyr_L,
#ਨਹੀਂਤਾਂ
  ssyr_U, ssyr_L,
#ਜੇ_ਅੰਤ
};

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*syr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  qsyr_thread_U, qsyr_thread_L,
#ਜਾਂਜੇ defined(DOUBLE)
  dsyr_thread_U, dsyr_thread_L,
#ਨਹੀਂਤਾਂ
  ssyr_thread_U, ssyr_thread_L,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  ਅੱਖਰ uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

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

  ਜੇ (lda  < MAX(1, n))   info =  7;
  ਜੇ (incx == 0)          info =  5;
  ਜੇ (n < 0)              info =  2;
  ਜੇ (uplo  < 0)          info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_UPLO Uplo, blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *a, blasint lda) {

  FLOAT *ਬਫਰ;
  ਪੂਰਨ uplo;
  blasint info;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 1;

    info = -1;

    ਜੇ (lda  < MAX(1, n))   info =  7;
    ਜੇ (incx == 0)          info =  5;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;

  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 1;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 0;

    info = -1;

    ਜੇ (lda  < MAX(1, n))   info =  7;
    ਜੇ (incx == 0)          info =  5;
    ਜੇ (n < 0)              info =  2;
    ਜੇ (uplo  < 0)          info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  ਜੇ (n == 0) ਵਾਪਸ;

  ਜੇ (alpha == ZERO) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();
#ਜੇ 1
  ਜੇ (incx == 1 && n < 100) {
    BLASLONG i;

    ਜੇ (uplo == 0) {
      ਲਈ (i = 0; i < n; i++){
        ਜੇ (x[i] != ZERO) {
          AXPYU_K(i + 1, 0, 0, alpha * x[i], x,     1, a, 1, ਨਿਰਾ, 0);
        }
        a += lda;
      }  
    } ਨਹੀਂਤਾਂ {
      ਲਈ (i = 0; i < n; i++){
        ਜੇ (x[i] != ZERO) {
          AXPYU_K(n - i, 0, 0, alpha * x[i], x + i, 1, a, 1, ਨਿਰਾ, 0);
        }
        a += 1 + lda;
      }
    }
    ਵਾਪਸ;
  } 
#ਜੇ_ਅੰਤ
  ਜੇ (incx < 0 ) x -= (n - 1) * incx;

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  nthreads = ਅੰਕ_cpu_avail(2);

  ਜੇ (n < 50) nthreads = 1;
  ਜੇ (nthreads > 2 && n < 500) nthreads = 2;

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  (syr[uplo])(n, alpha, x, incx, a, lda, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    (syr_thread[uplo])(n, alpha, x, incx, a, lda, ਬਫਰ, nthreads);

  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(1, n * n / 2 + n, n * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
