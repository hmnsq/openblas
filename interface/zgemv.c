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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#không_thì
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
static số_nguyên (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, số_nguyên) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#elif defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#không_thì
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(char *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  char trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  số_nguyên buffer_size;
#ifdef SMP
  số_nguyên nthreads;
#endif

  số_nguyên (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  nếu (trans == 'N')  i = 0;
  nếu (trans == 'T')  i = 1;
  nếu (trans == 'R')  i = 2;
  nếu (trans == 'C')  i = 3;
  nếu (trans == 'O')  i = 4;
  nếu (trans == 'U')  i = 5;
  nếu (trans == 'S')  i = 6;
  nếu (trans == 'D')  i = 7;

  nếu (incy == 0)      info = 11;
  nếu (incx == 0)      info = 8;
  nếu (lda < MAX(1,m)) info = 6;
  nếu (n < 0) 	      info = 3;
  nếu (m < 0) 	      info = 2;
  nếu (i < 0)          info = 1;

  trans = i;

  nếu (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#không_thì

void CNAME(enum CBLAS_ORDER order,
	   enum CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   void *VALPHA,
	   void  *va, blasint lda,
	   void  *vx, blasint incx,
	   void *VBETA,
	   void  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  số_nguyên trans, buffer_size;
  blasint info, t;
#ifdef SMP
  số_nguyên nthreads;
#endif

  số_nguyên (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  nếu (order == CblasColMajor) {
    nếu (TransA == CblasNoTrans)     trans = 0;
    nếu (TransA == CblasTrans)       trans = 1;
    nếu (TransA == CblasConjNoTrans) trans = 2;
    nếu (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    nếu (incy == 0)	  info = 11;
    nếu (incx == 0)	  info = 8;
    nếu (lda < MAX(1, m))  info = 6;
    nếu (n < 0)		  info = 3;
    nếu (m < 0)		  info = 2;
    nếu (trans < 0)        info = 1;

  }

  nếu (order == CblasRowMajor) {
    nếu (TransA == CblasNoTrans)     trans = 1;
    nếu (TransA == CblasTrans)       trans = 0;
    nếu (TransA == CblasConjNoTrans) trans = 3;
    nếu (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    nếu (incy == 0)	  info = 11;
    nếu (incx == 0)	  info = 8;
    nếu (lda < MAX(1, m))  info = 6;
    nếu (n < 0)		  info = 3;
    nếu (m < 0)		  info = 2;
    nếu (trans < 0)        info = 1;

  }

  nếu (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    trả_về;
  }

#endif

  /*  Quick return if possible. */

  nếu (m == 0 || n == 0) trả_về;

  lenx = n;
  leny = m;

  nếu (trans & 1) lenx = m;
  nếu (trans & 1) leny = n;

  nếu (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), NULL, 0, NULL, 0);

  nếu (alpha_r == ZERO && alpha_i == ZERO) trả_về;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  nếu (incx < 0) x -= (lenx - 1) * incx * 2;
  nếu (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / sizeof(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / sizeof(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  STACK_ALLOC(buffer_size, FLOAT, buffer);

#nếu defined(ARCH_X86_64) && defined(MAX_STACK_ALLOC) && MAX_STACK_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  nếu(trans && stack_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, sizeof(FLOAT) * buffer_size));
#endif

#ifdef SMP

#nếu defined(_WIN64) && defined(_M_ARM64)
  nếu (m*n > 25000000L) 
    nthreads = num_cpu_avail(4);
  không_thì
    nthreads = 1;
#không_thì
  nếu (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  không_thì
    nthreads = num_cpu_avail(2);
#endif  

  nếu (nthreads == 1) {
#endif  

    (gemv[(số_nguyên)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } không_thì {
    (gemv_thread[(số_nguyên)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  STACK_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  trả_về;
}
