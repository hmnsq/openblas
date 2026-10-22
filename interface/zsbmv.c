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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

/*
#ifdef SMP
#ifdef __64BIT__
#define SMPTEST 1
#endif
#endif
*/

#ifdef XDOUBLE
#define ERROR_NAME "XSBMV "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSBMV "
#وإلا
#define ERROR_NAME "CSBMV "
#endif

ساكن  صحيح (*sbmv[])(BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#ifdef XDOUBLE
  xsbmv_U, xsbmv_L,
#وإلاإذا defined(DOUBLE)
  zsbmv_U, zsbmv_L,
#وإلا
  csbmv_U, csbmv_L,
#endif
};

#ifdef SMPTEST
ساكن  صحيح (*sbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#ifdef XDOUBLE
  xsbmv_thread_U, xsbmv_thread_L,
#وإلاإذا defined(DOUBLE)
  zsbmv_thread_U, zsbmv_thread_L,
#وإلا
  csbmv_thread_U, csbmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
فراغ NAME(حرف *UPLO, blasint *N, blasint *K, FLOAT  *ALPHA, FLOAT *a, blasint *LDA,
            FLOAT  *b, blasint *INCX, FLOAT *BETA, FLOAT *c, blasint *INCY){

  حرف uplo_arg = *UPLO;
  blasint n	= *N;
  blasint k	= *K;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;
  FLOAT beta_r	= BETA[0];
  FLOAT beta_i	= BETA[1];
  blasint incy	= *INCY;

  blasint info;
  صحيح uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  صحيح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  𐤀𐤌 (uplo_arg  == 'U') uplo  = 0;
  𐤀𐤌 (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  𐤀𐤌 (incy == 0)          info = 11;
  𐤀𐤌 (incx == 0)          info =  8;
  𐤀𐤌 (lda  < k + 1)       info =  6;
  𐤀𐤌 (k < 0)              info =  3;
  𐤀𐤌 (n < 0)              info =  2;
  𐤀𐤌 (uplo  < 0)          info =  1;

  𐤀𐤌 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  𐤀𐤌 (n == 0) أرجع;

  𐤀𐤌 ((beta_r != ONE) || (beta_i != ZERO)) SCAL_K(n, 0, 0, beta_r, beta_i, c, blasabs(incy), عدم, 0, عدم, 0);

  𐤀𐤌 ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  𐤀𐤌 (incx < 0 ) b -= (n - 1) * incx * COMPSIZE;
  𐤀𐤌 (incy < 0 ) c -= (n - 1) * incy * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  𐤀𐤌 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#ifdef SMPTEST
  nthreads = num_cpu_avail(2);

  𐤀𐤌 (nthreads == 1) {
#endif

  (sbmv[uplo])(n, k, alpha_r, alpha_i, a, lda, b, incx, c, incy, buffer);

#ifdef SMPTEST
  } وإلا {

    (sbmv_thread[uplo])(n, k, ALPHA, a, lda, b, incx, c, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  أرجع;
}
