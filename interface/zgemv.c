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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, പൂർണ്ണസംഖ്യ) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#അല്ലെങ്കിലെങ്കിൽ defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#അല്ലെങ്കിൽ
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  അക്ഷരം trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  പൂർണ്ണസംഖ്യ buffer_size;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  എങ്കിൽ (trans == 'N')  i = 0;
  എങ്കിൽ (trans == 'T')  i = 1;
  എങ്കിൽ (trans == 'R')  i = 2;
  എങ്കിൽ (trans == 'C')  i = 3;
  എങ്കിൽ (trans == 'O')  i = 4;
  എങ്കിൽ (trans == 'U')  i = 5;
  എങ്കിൽ (trans == 'S')  i = 6;
  എങ്കിൽ (trans == 'D')  i = 7;

  എങ്കിൽ (incy == 0)      info = 11;
  എങ്കിൽ (incx == 0)      info = 8;
  എങ്കിൽ (lda < MAX(1,m)) info = 6;
  എങ്കിൽ (n < 0) 	      info = 3;
  എങ്കിൽ (m < 0) 	      info = 2;
  എങ്കിൽ (i < 0)          info = 1;

  trans = i;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order,
	   എണ്ണൽ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   ശൂന്യം *VALPHA,
	   ശൂന്യം  *va, blasint lda,
	   ശൂന്യം  *vx, blasint incx,
	   ശൂന്യം *VBETA,
	   ശൂന്യം  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  പൂർണ്ണസംഖ്യ trans, buffer_size;
  blasint info, t;
#ifdef SMP
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 0;
    എങ്കിൽ (TransA == CblasTrans)       trans = 1;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 2;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    എങ്കിൽ (incy == 0)	  info = 11;
    എങ്കിൽ (incx == 0)	  info = 8;
    എങ്കിൽ (lda < MAX(1, m))  info = 6;
    എങ്കിൽ (n < 0)		  info = 3;
    എങ്കിൽ (m < 0)		  info = 2;
    എങ്കിൽ (trans < 0)        info = 1;

  }

  എങ്കിൽ (order == CblasRowMajor) {
    എങ്കിൽ (TransA == CblasNoTrans)     trans = 1;
    എങ്കിൽ (TransA == CblasTrans)       trans = 0;
    എങ്കിൽ (TransA == CblasConjNoTrans) trans = 3;
    എങ്കിൽ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    എങ്കിൽ (incy == 0)	  info = 11;
    എങ്കിൽ (incx == 0)	  info = 8;
    എങ്കിൽ (lda < MAX(1, m))  info = 6;
    എങ്കിൽ (n < 0)		  info = 3;
    എങ്കിൽ (m < 0)		  info = 2;
    എങ്കിൽ (trans < 0)        info = 1;

  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#endif

  /*  Quick return if possible. */

  എങ്കിൽ (m == 0 || n == 0) തിരികെ;

  lenx = n;
  leny = m;

  എങ്കിൽ (trans & 1) lenx = m;
  എങ്കിൽ (trans & 1) leny = n;

  എങ്കിൽ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ശൂന്യ, 0, ശൂന്യ, 0);

  എങ്കിൽ (alpha_r == ZERO && alpha_i == ZERO) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  എങ്കിൽ (incx < 0) x -= (lenx - 1) * incx * 2;
  എങ്കിൽ (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / വലുപ്പം(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / വലുപ്പം(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  അട്ടി_ALLOC(buffer_size, FLOAT, buffer);

#എങ്കിൽ defined(ARCH_X86_64) && defined(MAX_അട്ടി_ALLOC) && MAX_അട്ടി_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  എങ്കിൽ(trans && അട്ടി_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, വലുപ്പം(FLOAT) * buffer_size));
#endif

#ifdef SMP

#എങ്കിൽ defined(_WIN64) && defined(_M_ARM64)
  എങ്കിൽ (m*n > 25000000L) 
    nthreads = അക്കം_cpu_avail(4);
  അല്ലെങ്കിൽ
    nthreads = 1;
#അല്ലെങ്കിൽ
  എങ്കിൽ (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  അല്ലെങ്കിൽ
    nthreads = അക്കം_cpu_avail(2);
#endif  

  എങ്കിൽ (nthreads == 1) {
#endif  

    (gemv[(പൂർണ്ണസംഖ്യ)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } അല്ലെങ്കിൽ {
    (gemv_thread[(പൂർണ്ണസംഖ്യ)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  അട്ടി_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  തിരികെ;
}
