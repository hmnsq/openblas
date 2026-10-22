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
#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGEMV "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGEMV "
#وإلا
#وصف_ڪريو ERROR_NAME "CGEMV "
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
ساكن صحيح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحيح) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#وإلاإذا defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#وإلا
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#جيڪڏهن_پورو
};
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  حرف trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *بفر;
  صحيح بفر_سائيز;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  صحيح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  إذا (trans == 'N')  i = 0;
  إذا (trans == 'T')  i = 1;
  إذا (trans == 'R')  i = 2;
  إذا (trans == 'C')  i = 3;
  إذا (trans == 'O')  i = 4;
  إذا (trans == 'U')  i = 5;
  إذا (trans == 'S')  i = 6;
  إذا (trans == 'D')  i = 7;

  إذا (incy == 0)      info = 11;
  إذا (incx == 0)      info = 8;
  إذا (lda < MAX(1,m)) info = 6;
  إذا (n < 0) 	      info = 3;
  إذا (m < 0) 	      info = 2;
  إذا (i < 0)          info = 1;

  trans = i;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب,
	   تعداد CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   فراغ *VALPHA,
	   فراغ  *va, blasint lda,
	   فراغ  *vx, blasint incx,
	   فراغ *VBETA,
	   فراغ  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *بفر;
  blasint    lenx, leny;
  صحيح trans, بفر_سائيز;
  blasint info, t;
#جيڪڏهن_وصف_ٿيل SMP
  صحيح nthreads;
#جيڪڏهن_پورو

  صحيح (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 2;
    إذا (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    إذا (incy == 0)	  info = 11;
    إذا (incx == 0)	  info = 8;
    إذا (lda < MAX(1, m))  info = 6;
    إذا (n < 0)		  info = 3;
    إذا (m < 0)		  info = 2;
    إذا (trans < 0)        info = 1;

  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 3;
    إذا (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    إذا (incy == 0)	  info = 11;
    إذا (incx == 0)	  info = 8;
    إذا (lda < MAX(1, m))  info = 6;
    إذا (n < 0)		  info = 3;
    إذا (m < 0)		  info = 2;
    إذا (trans < 0)        info = 1;

  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_پورو

  /*  Quick return if possible. */

  إذا (m == 0 || n == 0) أرجع;

  lenx = n;
  leny = m;

  إذا (trans & 1) lenx = m;
  إذا (trans & 1) leny = n;

  إذا (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), عدم, 0, عدم, 0);

  إذا (alpha_r == ZERO && alpha_i == ZERO) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (lenx - 1) * incx * 2;
  إذا (incy < 0) y -= (leny - 1) * incy * 2;

  بفر_سائيز = 2 * (m + n) + 128 / حجم(FLOAT);
#جيڪڏهن_وصف_ٿيل WINDOWS_ABI
  بفر_سائيز += 160 / حجم(FLOAT) ;
#جيڪڏهن_پورو
  // for alignment
  بفر_سائيز = (بفر_سائيز + 3) & ~3;
  ڍير_ALLOC(بفر_سائيز, FLOAT, بفر);

#إذا defined(ARCH_X86_64) && defined(MAX_ڍير_ALLOC) && MAX_ڍير_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  إذا(trans && ڍير_alloc_سائيز)
    memset(بفر, 0, MIN(بفر_سائيز_2, حجم(FLOAT) * بفر_سائيز));
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP

#إذا defined(_WIN64) && defined(_M_ARM64)
  إذا (m*n > 25000000L) 
    nthreads = انگ_cpu_avail(4);
  وإلا
    nthreads = 1;
#وإلا
  إذا (1L * m * n < 1024L * GEMM_MULTITHREAD_حد)
    nthreads = 1;
  وإلا
    nthreads = انگ_cpu_avail(2);
#جيڪڏهن_پورو  

  إذا (nthreads == 1) {
#جيڪڏهن_پورو  

    (gemv[(صحيح)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, بفر);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {
    (gemv_thread[(صحيح)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, بفر, nthreads);
  }
#جيڪڏهن_پورو


  ڍير_FREE(بفر);

  FUNCTION_PROFILE_پڇاڙي(4, m * n + m + n,  2 * m * n);

  IDEBUG_پڇاڙي;

  أرجع;
}
