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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ପୂର୍ଣ୍ଣସଂଖ୍ୟା) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#ନଚେତଯଦି defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#ଅନ୍ୟଥା
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(ଅକ୍ଷର *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  ଅକ୍ଷର trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା buffer_size;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ଯଦି (trans == 'N')  i = 0;
  ଯଦି (trans == 'T')  i = 1;
  ଯଦି (trans == 'R')  i = 2;
  ଯଦି (trans == 'C')  i = 3;
  ଯଦି (trans == 'O')  i = 4;
  ଯଦି (trans == 'U')  i = 5;
  ଯଦି (trans == 'S')  i = 6;
  ଯଦି (trans == 'D')  i = 7;

  ଯଦି (incy == 0)      info = 11;
  ଯଦି (incx == 0)      info = 8;
  ଯଦି (lda < MAX(1,m)) info = 6;
  ଯଦି (n < 0) 	      info = 3;
  ଯଦି (m < 0) 	      info = 2;
  ଯଦି (i < 0)          info = 1;

  trans = i;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order,
	   ଗଣନା CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   ଶୂନ୍ୟ *VALPHA,
	   ଶୂନ୍ୟ  *va, blasint lda,
	   ଶୂନ୍ୟ  *vx, blasint incx,
	   ଶୂନ୍ୟ *VBETA,
	   ଶୂନ୍ୟ  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans, buffer_size;
  blasint info, t;
#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା nthreads;
#endif

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ଯଦି (order == CblasColMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 0;
    ଯଦି (TransA == CblasTrans)       trans = 1;
    ଯଦି (TransA == CblasConjNoTrans) trans = 2;
    ଯଦି (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ଯଦି (incy == 0)	  info = 11;
    ଯଦି (incx == 0)	  info = 8;
    ଯଦି (lda < MAX(1, m))  info = 6;
    ଯଦି (n < 0)		  info = 3;
    ଯଦି (m < 0)		  info = 2;
    ଯଦି (trans < 0)        info = 1;

  }

  ଯଦି (order == CblasRowMajor) {
    ଯଦି (TransA == CblasNoTrans)     trans = 1;
    ଯଦି (TransA == CblasTrans)       trans = 0;
    ଯଦି (TransA == CblasConjNoTrans) trans = 3;
    ଯଦି (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    ଯଦି (incy == 0)	  info = 11;
    ଯଦି (incx == 0)	  info = 8;
    ଯଦି (lda < MAX(1, m))  info = 6;
    ଯଦି (n < 0)		  info = 3;
    ଯଦି (m < 0)		  info = 2;
    ଯଦି (trans < 0)        info = 1;

  }

  ଯଦି (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
    ଫେରନ୍ତୁ;
  }

#endif

  /*  Quick return if possible. */

  ଯଦି (m == 0 || n == 0) ଫେରନ୍ତୁ;

  lenx = n;
  leny = m;

  ଯଦି (trans & 1) lenx = m;
  ଯଦି (trans & 1) leny = n;

  ଯଦି (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ଖାଲି, 0, ଖାଲି, 0);

  ଯଦି (alpha_r == ZERO && alpha_i == ZERO) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0) x -= (lenx - 1) * incx * 2;
  ଯଦି (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / ଆକାର(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / ଆକାର(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ଗଦା_ALLOC(buffer_size, FLOAT, buffer);

#ଯଦି defined(ARCH_X86_64) && defined(MAX_ଗଦା_ALLOC) && MAX_ଗଦା_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  ଯଦି(trans && ଗଦା_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, ଆକାର(FLOAT) * buffer_size));
#endif

#ifdef SMP

#ଯଦି defined(_WIN64) && defined(_M_ARM64)
  ଯଦି (m*n > 25000000L) 
    nthreads = ଅଙ୍କ_cpu_avail(4);
  ଅନ୍ୟଥା
    nthreads = 1;
#ଅନ୍ୟଥା
  ଯଦି (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  ଅନ୍ୟଥା
    nthreads = ଅଙ୍କ_cpu_avail(2);
#endif  

  ଯଦି (nthreads == 1) {
#endif  

    (gemv[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ଅନ୍ୟଥା {
    (gemv_thread[(ପୂର୍ଣ୍ଣସଂଖ୍ୟା)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  ଗଦା_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
