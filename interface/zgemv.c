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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#નહીંતર
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
સ્થિર પૂર્ણાંક (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#અથવાજો defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#નહીંતર
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  અક્ષર trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  પૂર્ણાંક buffer_size;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  જો (trans == 'N')  i = 0;
  જો (trans == 'T')  i = 1;
  જો (trans == 'R')  i = 2;
  જો (trans == 'C')  i = 3;
  જો (trans == 'O')  i = 4;
  જો (trans == 'U')  i = 5;
  જો (trans == 'S')  i = 6;
  જો (trans == 'D')  i = 7;

  જો (incy == 0)      info = 11;
  જો (incx == 0)      info = 8;
  જો (lda < MAX(1,m)) info = 6;
  જો (n < 0) 	      info = 3;
  જો (m < 0) 	      info = 2;
  જો (i < 0)          info = 1;

  trans = i;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   ગણના CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   શૂન્ય *VALPHA,
	   શૂન્ય  *va, blasint lda,
	   શૂન્ય  *vx, blasint incx,
	   શૂન્ય *VBETA,
	   શૂન્ય  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  પૂર્ણાંક trans, buffer_size;
  blasint info, t;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  જો (order == CblasColMajor) {
    જો (TransA == CblasNoTrans)     trans = 0;
    જો (TransA == CblasTrans)       trans = 1;
    જો (TransA == CblasConjNoTrans) trans = 2;
    જો (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    જો (incy == 0)	  info = 11;
    જો (incx == 0)	  info = 8;
    જો (lda < MAX(1, m))  info = 6;
    જો (n < 0)		  info = 3;
    જો (m < 0)		  info = 2;
    જો (trans < 0)        info = 1;

  }

  જો (order == CblasRowMajor) {
    જો (TransA == CblasNoTrans)     trans = 1;
    જો (TransA == CblasTrans)       trans = 0;
    જો (TransA == CblasConjNoTrans) trans = 3;
    જો (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    જો (incy == 0)	  info = 11;
    જો (incx == 0)	  info = 8;
    જો (lda < MAX(1, m))  info = 6;
    જો (n < 0)		  info = 3;
    જો (m < 0)		  info = 2;
    જો (trans < 0)        info = 1;

  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  /*  Quick return if possible. */

  જો (m == 0 || n == 0) પાછા;

  lenx = n;
  leny = m;

  જો (trans & 1) lenx = m;
  જો (trans & 1) leny = n;

  જો (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ખાલી, 0, ખાલી, 0);

  જો (alpha_r == ZERO && alpha_i == ZERO) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (lenx - 1) * incx * 2;
  જો (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / કદ(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / કદ(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ઢગલો_ALLOC(buffer_size, FLOAT, buffer);

#જો defined(ARCH_X86_64) && defined(MAX_ઢગલો_ALLOC) && MAX_ઢગલો_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  જો(trans && ઢગલો_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, કદ(FLOAT) * buffer_size));
#endif

#ifdef SMP

#જો defined(_WIN64) && defined(_M_ARM64)
  જો (m*n > 25000000L) 
    nthreads = અંક_cpu_avail(4);
  નહીંતર
    nthreads = 1;
#નહીંતર
  જો (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  નહીંતર
    nthreads = અંક_cpu_avail(2);
#endif  

  જો (nthreads == 1) {
#endif  

    (gemv[(પૂર્ણાંક)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } નહીંતર {
    (gemv_thread[(પૂર્ણાંક)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  ઢગલો_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  પાછા;
}
