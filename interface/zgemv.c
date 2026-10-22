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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಙ್ಕ) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#ನೋಚೇದ್ಯದಿ defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#ಅನ್ಯಥಾ
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  ವರ್ಣ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  ಪೂರ್ಣಾಙ್ಕ buffer_size;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ಯದಿ (trans == 'N')  i = 0;
  ಯದಿ (trans == 'T')  i = 1;
  ಯದಿ (trans == 'R')  i = 2;
  ಯದಿ (trans == 'C')  i = 3;
  ಯದಿ (trans == 'O')  i = 4;
  ಯದಿ (trans == 'U')  i = 5;
  ಯದಿ (trans == 'S')  i = 6;
  ಯದಿ (trans == 'D')  i = 7;

  ಯದಿ (incy == 0)      info = 11;
  ಯದಿ (incx == 0)      info = 8;
  ಯದಿ (lda < MAX(1,m)) info = 6;
  ಯದಿ (n < 0) 	      info = 3;
  ಯದಿ (m < 0) 	      info = 2;
  ಯದಿ (i < 0)          info = 1;

  trans = i;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order,
	   ಗಣನಾ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   ಶೂನ್ಯ *VALPHA,
	   ಶೂನ್ಯ  *va, blasint lda,
	   ಶೂನ್ಯ  *vx, blasint incx,
	   ಶೂನ್ಯ *VBETA,
	   ಶೂನ್ಯ  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  ಪೂರ್ಣಾಙ್ಕ trans, buffer_size;
  blasint info, t;
#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

  ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 0;
    ಯದಿ (TransA == CblasTrans)       trans = 1;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 2;
    ಯದಿ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ಯದಿ (incy == 0)	  info = 11;
    ಯದಿ (incx == 0)	  info = 8;
    ಯದಿ (lda < MAX(1, m))  info = 6;
    ಯದಿ (n < 0)		  info = 3;
    ಯದಿ (m < 0)		  info = 2;
    ಯದಿ (trans < 0)        info = 1;

  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (TransA == CblasNoTrans)     trans = 1;
    ಯದಿ (TransA == CblasTrans)       trans = 0;
    ಯದಿ (TransA == CblasConjNoTrans) trans = 3;
    ಯದಿ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    ಯದಿ (incy == 0)	  info = 11;
    ಯದಿ (incx == 0)	  info = 8;
    ಯದಿ (lda < MAX(1, m))  info = 6;
    ಯದಿ (n < 0)		  info = 3;
    ಯದಿ (m < 0)		  info = 2;
    ಯದಿ (trans < 0)        info = 1;

  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#endif

  /*  Quick return if possible. */

  ಯದಿ (m == 0 || n == 0) ಪ್ರತ್ಯಾಗಮ;

  lenx = n;
  leny = m;

  ಯದಿ (trans & 1) lenx = m;
  ಯದಿ (trans & 1) leny = n;

  ಯದಿ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);

  ಯದಿ (alpha_r == ZERO && alpha_i == ZERO) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (lenx - 1) * incx * 2;
  ಯದಿ (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / ಪರಿಮಾಣ(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / ಪರಿಮಾಣ(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  ರಾಶಿ_ALLOC(buffer_size, FLOAT, buffer);

#ಯದಿ defined(ARCH_X86_64) && defined(MAX_ರಾಶಿ_ALLOC) && MAX_ರಾಶಿ_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  ಯದಿ(trans && ರಾಶಿ_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, ಪರಿಮಾಣ(FLOAT) * buffer_size));
#endif

#ifdef SMP

#ಯದಿ defined(_WIN64) && defined(_M_ARM64)
  ಯದಿ (m*n > 25000000L) 
    nthreads = ಅಂಕೆ_cpu_avail(4);
  ಅನ್ಯಥಾ
    nthreads = 1;
#ಅನ್ಯಥಾ
  ಯದಿ (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  ಅನ್ಯಥಾ
    nthreads = ಅಂಕೆ_cpu_avail(2);
#endif  

  ಯದಿ (nthreads == 1) {
#endif  

    (gemv[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } ಅನ್ಯಥಾ {
    (gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  ರಾಶಿ_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
