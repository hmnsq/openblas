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
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#নতুবা_যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#অন্যথা
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
স্থির পূর্ণসংখ্যা (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, পূর্ণসংখ্যা) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#নতুবা_যদি defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#অন্যথা
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(অক্ষর *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  অক্ষর trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *বাফার;
  পূর্ণসংখ্যা বাফার_আকার;
#ifdef SMP
  পূর্ণসংখ্যা nthreads;
#endif

  পূর্ণসংখ্যা (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  যদি (trans == 'N')  i = 0;
  যদি (trans == 'T')  i = 1;
  যদি (trans == 'R')  i = 2;
  যদি (trans == 'C')  i = 3;
  যদি (trans == 'O')  i = 4;
  যদি (trans == 'U')  i = 5;
  যদি (trans == 'S')  i = 6;
  যদি (trans == 'D')  i = 7;

  যদি (incy == 0)      info = 11;
  যদি (incx == 0)      info = 8;
  যদি (lda < MAX(1,m)) info = 6;
  যদি (n < 0) 	      info = 3;
  যদি (m < 0) 	      info = 2;
  যদি (i < 0)          info = 1;

  trans = i;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order,
	   গণনা CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   শূন্য *VALPHA,
	   শূন্য  *va, blasint lda,
	   শূন্য  *vx, blasint incx,
	   শূন্য *VBETA,
	   শূন্য  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *বাফার;
  blasint    lenx, leny;
  পূর্ণসংখ্যা trans, বাফার_আকার;
  blasint info, t;
#ifdef SMP
  পূর্ণসংখ্যা nthreads;
#endif

  পূর্ণসংখ্যা (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  যদি (order == CblasColMajor) {
    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 2;
    যদি (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    যদি (incy == 0)	  info = 11;
    যদি (incx == 0)	  info = 8;
    যদি (lda < MAX(1, m))  info = 6;
    যদি (n < 0)		  info = 3;
    যদি (m < 0)		  info = 2;
    যদি (trans < 0)        info = 1;

  }

  যদি (order == CblasRowMajor) {
    যদি (TransA == CblasNoTrans)     trans = 1;
    যদি (TransA == CblasTrans)       trans = 0;
    যদি (TransA == CblasConjNoTrans) trans = 3;
    যদি (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    যদি (incy == 0)	  info = 11;
    যদি (incx == 0)	  info = 8;
    যদি (lda < MAX(1, m))  info = 6;
    যদি (n < 0)		  info = 3;
    যদি (m < 0)		  info = 2;
    যদি (trans < 0)        info = 1;

  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, আকার(ERROR_NAME));
    ফেরত;
  }

#endif

  /*  Quick return if possible. */

  যদি (m == 0 || n == 0) ফেরত;

  lenx = n;
  leny = m;

  যদি (trans & 1) lenx = m;
  যদি (trans & 1) leny = n;

  যদি (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ফাঁকা, 0, ফাঁকা, 0);

  যদি (alpha_r == ZERO && alpha_i == ZERO) ফেরত;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (lenx - 1) * incx * 2;
  যদি (incy < 0) y -= (leny - 1) * incy * 2;

  বাফার_আকার = 2 * (m + n) + 128 / আকার(FLOAT);
#ifdef WINDOWS_ABI
  বাফার_আকার += 160 / আকার(FLOAT) ;
#endif
  // for alignment
  বাফার_আকার = (বাফার_আকার + 3) & ~3;
  স্তূপ_ALLOC(বাফার_আকার, FLOAT, বাফার);

#যদি defined(ARCH_X86_64) && defined(MAX_স্তূপ_ALLOC) && MAX_স্তূপ_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  যদি(trans && স্তূপ_alloc_আকার)
    memset(বাফার, 0, MIN(বাফার_আকার_2, আকার(FLOAT) * বাফার_আকার));
#endif

#ifdef SMP

#যদি defined(_WIN64) && defined(_M_ARM64)
  যদি (m*n > 25000000L) 
    nthreads = অঙ্ক_cpu_avail(4);
  অন্যথা
    nthreads = 1;
#অন্যথা
  যদি (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  অন্যথা
    nthreads = অঙ্ক_cpu_avail(2);
#endif  

  যদি (nthreads == 1) {
#endif  

    (gemv[(পূর্ণসংখ্যা)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, বাফার);

#ifdef SMP
  } অন্যথা {
    (gemv_thread[(পূর্ণসংখ্যা)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, বাফার, nthreads);
  }
#endif


  স্তূপ_FREE(বাফার);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ফেরত;
}
