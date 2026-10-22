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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGEMV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#अञ्ञथा
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
स्थिर पूर्णाङ्क (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#नोचेद्यदि defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#अञ्ञथा
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  वर्ण trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  पूर्णाङ्क buffer_size;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  सचे (trans == 'N')  i = 0;
  सचे (trans == 'T')  i = 1;
  सचे (trans == 'R')  i = 2;
  सचे (trans == 'C')  i = 3;
  सचे (trans == 'O')  i = 4;
  सचे (trans == 'U')  i = 5;
  सचे (trans == 'S')  i = 6;
  सचे (trans == 'D')  i = 7;

  सचे (incy == 0)      info = 11;
  सचे (incx == 0)      info = 8;
  सचे (lda < MAX(1,m)) info = 6;
  सचे (n < 0) 	      info = 3;
  सचे (m < 0) 	      info = 2;
  सचे (i < 0)          info = 1;

  trans = i;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   शून्य *VALPHA,
	   शून्य  *va, blasint lda,
	   शून्य  *vx, blasint incx,
	   शून्य *VBETA,
	   शून्य  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *buffer;
  blasint    lenx, leny;
  पूर्णाङ्क trans, buffer_size;
  blasint info, t;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  सचे (order == CblasColMajor) {
    सचे (TransA == CblasNoTrans)     trans = 0;
    सचे (TransA == CblasTrans)       trans = 1;
    सचे (TransA == CblasConjNoTrans) trans = 2;
    सचे (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    सचे (incy == 0)	  info = 11;
    सचे (incx == 0)	  info = 8;
    सचे (lda < MAX(1, m))  info = 6;
    सचे (n < 0)		  info = 3;
    सचे (m < 0)		  info = 2;
    सचे (trans < 0)        info = 1;

  }

  सचे (order == CblasRowMajor) {
    सचे (TransA == CblasNoTrans)     trans = 1;
    सचे (TransA == CblasTrans)       trans = 0;
    सचे (TransA == CblasConjNoTrans) trans = 3;
    सचे (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    सचे (incy == 0)	  info = 11;
    सचे (incx == 0)	  info = 8;
    सचे (lda < MAX(1, m))  info = 6;
    सचे (n < 0)		  info = 3;
    सचे (m < 0)		  info = 2;
    सचे (trans < 0)        info = 1;

  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  /*  Quick return if possible. */

  सचे (m == 0 || n == 0) पटिनिवत्त;

  lenx = n;
  leny = m;

  सचे (trans & 1) lenx = m;
  सचे (trans & 1) leny = n;

  सचे (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  सचे (alpha_r == ZERO && alpha_i == ZERO) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0) x -= (lenx - 1) * incx * 2;
  सचे (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / परिमाण(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / परिमाण(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  राशिः_ALLOC(buffer_size, FLOAT, buffer);

#सचे defined(ARCH_X86_64) && defined(MAX_राशिः_ALLOC) && MAX_राशिः_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  सचे(trans && राशिः_alloc_size)
    memset(buffer, 0, MIN(BUFFER_SIZE, परिमाण(FLOAT) * buffer_size));
#endif

#ifdef SMP

#सचे defined(_WIN64) && defined(_M_ARM64)
  सचे (m*n > 25000000L) 
    nthreads = अङ्कः_cpu_avail(4);
  अञ्ञथा
    nthreads = 1;
#अञ्ञथा
  सचे (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  अञ्ञथा
    nthreads = अङ्कः_cpu_avail(2);
#endif  

  सचे (nthreads == 1) {
#endif  

    (gemv[(पूर्णाङ्क)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अञ्ञथा {
    (gemv_thread[(पूर्णाङ्क)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  राशिः_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  पटिनिवत्त;
}
