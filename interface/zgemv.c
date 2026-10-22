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
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMV "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMV "
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ਪੂਰਨ) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#ਜਾਂਜੇ defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#ਨਹੀਂਤਾਂ
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#ਜੇ_ਅੰਤ
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  ਅੱਖਰ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *ਬਫਰ;
  ਪੂਰਨ ਬਫਰ_ਮਾਪ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
		FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
		  GEMV_N, GEMV_T, GEMV_R, GEMV_C,
		  GEMV_O, GEMV_U, GEMV_S, GEMV_D,
		};

  blasint    info;
  blasint    lenx, leny;
  blasint    i;

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

  FLOAT beta_r  = *(BETA + 0);
  FLOAT beta_i  = *(BETA + 1);

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i    = -1;

  ਜੇ (trans == 'N')  i = 0;
  ਜੇ (trans == 'T')  i = 1;
  ਜੇ (trans == 'R')  i = 2;
  ਜੇ (trans == 'C')  i = 3;
  ਜੇ (trans == 'O')  i = 4;
  ਜੇ (trans == 'U')  i = 5;
  ਜੇ (trans == 'S')  i = 6;
  ਜੇ (trans == 'D')  i = 7;

  ਜੇ (incy == 0)      info = 11;
  ਜੇ (incx == 0)      info = 8;
  ਜੇ (lda < MAX(1,m)) info = 6;
  ਜੇ (n < 0) 	      info = 3;
  ਜੇ (m < 0) 	      info = 2;
  ਜੇ (i < 0)          info = 1;

  trans = i;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   ਖਾਲੀ *VALPHA,
	   ਖਾਲੀ  *va, blasint lda,
	   ਖਾਲੀ  *vx, blasint incx,
	   ਖਾਲੀ *VBETA,
	   ਖਾਲੀ  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *ਬਫਰ;
  blasint    lenx, leny;
  ਪੂਰਨ trans, ਬਫਰ_ਮਾਪ;
  blasint info, t;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

  ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
	    FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
	      GEMV_N, GEMV_T, GEMV_R, GEMV_C,
	      GEMV_O, GEMV_U, GEMV_S, GEMV_D,
	    };

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

  FLOAT beta_r  = *(BETA + 0);
  FLOAT beta_i  = *(BETA + 1);

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 0;
    ਜੇ (TransA == CblasTrans)       trans = 1;
    ਜੇ (TransA == CblasConjNoTrans) trans = 2;
    ਜੇ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ਜੇ (incy == 0)	  info = 11;
    ਜੇ (incx == 0)	  info = 8;
    ਜੇ (lda < MAX(1, m))  info = 6;
    ਜੇ (n < 0)		  info = 3;
    ਜੇ (m < 0)		  info = 2;
    ਜੇ (trans < 0)        info = 1;

  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (TransA == CblasNoTrans)     trans = 1;
    ਜੇ (TransA == CblasTrans)       trans = 0;
    ਜੇ (TransA == CblasConjNoTrans) trans = 3;
    ਜੇ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    ਜੇ (incy == 0)	  info = 11;
    ਜੇ (incx == 0)	  info = 8;
    ਜੇ (lda < MAX(1, m))  info = 6;
    ਜੇ (n < 0)		  info = 3;
    ਜੇ (m < 0)		  info = 2;
    ਜੇ (trans < 0)        info = 1;

  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਜੇ_ਅੰਤ

  /*  Quick return if possible. */

  ਜੇ (m == 0 || n == 0) ਵਾਪਸ;

  lenx = n;
  leny = m;

  ਜੇ (trans & 1) lenx = m;
  ਜੇ (trans & 1) leny = n;

  ਜੇ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);

  ਜੇ (alpha_r == ZERO && alpha_i == ZERO) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (lenx - 1) * incx * 2;
  ਜੇ (incy < 0) y -= (leny - 1) * incy * 2;

  ਬਫਰ_ਮਾਪ = 2 * (m + n) + 128 / ਆਕਾਰ(FLOAT);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ WINDOWS_ABI
  ਬਫਰ_ਮਾਪ += 160 / ਆਕਾਰ(FLOAT) ;
#ਜੇ_ਅੰਤ
  // for alignment
  ਬਫਰ_ਮਾਪ = (ਬਫਰ_ਮਾਪ + 3) & ~3;
  ਢੇਰ_ALLOC(ਬਫਰ_ਮਾਪ, FLOAT, ਬਫਰ);

#ਜੇ defined(ARCH_X86_64) && defined(MAX_ਢੇਰ_ALLOC) && MAX_ਢੇਰ_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  ਜੇ(trans && ਢੇਰ_alloc_ਮਾਪ)
    memset(ਬਫਰ, 0, MIN(ਬਫਰ_ਮਾਪ_2, ਆਕਾਰ(FLOAT) * ਬਫਰ_ਮਾਪ));
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP

#ਜੇ defined(_WIN64) && defined(_M_ARM64)
  ਜੇ (m*n > 25000000L) 
    nthreads = ਅੰਕ_cpu_avail(4);
  ਨਹੀਂਤਾਂ
    nthreads = 1;
#ਨਹੀਂਤਾਂ
  ਜੇ (1L * m * n < 1024L * GEMM_MULTITHREAD_ਹੱਦ)
    nthreads = 1;
  ਨਹੀਂਤਾਂ
    nthreads = ਅੰਕ_cpu_avail(2);
#ਜੇ_ਅੰਤ  

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ  

    (gemv[(ਪੂਰਨ)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, ਬਫਰ);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {
    (gemv_thread[(ਪੂਰਨ)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, ਬਫਰ, nthreads);
  }
#ਜੇ_ਅੰਤ


  ਢੇਰ_FREE(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(4, m * n + m + n,  2 * m * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
