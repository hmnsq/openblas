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
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZGEMV "
#अन्यथा
#define ERROR_NAME "CGEMV "
#endif

#ifdef SMP
स्थिर पूर्णांक (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#ifdef XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#किंवाजर defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#अन्यथा
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *buffer;
  पूर्णांक buffer_size;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  जर (trans == 'N')  i = 0;
  जर (trans == 'T')  i = 1;
  जर (trans == 'R')  i = 2;
  जर (trans == 'C')  i = 3;
  जर (trans == 'O')  i = 4;
  जर (trans == 'U')  i = 5;
  जर (trans == 'S')  i = 6;
  जर (trans == 'D')  i = 7;

  जर (incy == 0)      info = 11;
  जर (incx == 0)      info = 8;
  जर (lda < MAX(1,m)) info = 6;
  जर (n < 0) 	      info = 3;
  जर (m < 0) 	      info = 2;
  जर (i < 0)          info = 1;

  trans = i;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

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
  पूर्णांक trans, buffer_size;
  blasint info, t;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  जर (order == CblasColMajor) {
    जर (TransA == CblasNoTrans)     trans = 0;
    जर (TransA == CblasTrans)       trans = 1;
    जर (TransA == CblasConjNoTrans) trans = 2;
    जर (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    जर (incy == 0)	  info = 11;
    जर (incx == 0)	  info = 8;
    जर (lda < MAX(1, m))  info = 6;
    जर (n < 0)		  info = 3;
    जर (m < 0)		  info = 2;
    जर (trans < 0)        info = 1;

  }

  जर (order == CblasRowMajor) {
    जर (TransA == CblasNoTrans)     trans = 1;
    जर (TransA == CblasTrans)       trans = 0;
    जर (TransA == CblasConjNoTrans) trans = 3;
    जर (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    जर (incy == 0)	  info = 11;
    जर (incx == 0)	  info = 8;
    जर (lda < MAX(1, m))  info = 6;
    जर (n < 0)		  info = 3;
    जर (m < 0)		  info = 2;
    जर (trans < 0)        info = 1;

  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#endif

  /*  Quick return if possible. */

  जर (m == 0 || n == 0) परत;

  lenx = n;
  leny = m;

  जर (trans & 1) lenx = m;
  जर (trans & 1) leny = n;

  जर (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  जर (alpha_r == ZERO && alpha_i == ZERO) परत;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जर (incx < 0) x -= (lenx - 1) * incx * 2;
  जर (incy < 0) y -= (leny - 1) * incy * 2;

  buffer_size = 2 * (m + n) + 128 / आकार(FLOAT);
#ifdef WINDOWS_ABI
  buffer_size += 160 / आकार(FLOAT) ;
#endif
  // for alignment
  buffer_size = (buffer_size + 3) & ~3;
  चळत_ALLOC(buffer_size, FLOAT, buffer);

#जर defined(ARCH_X86_64) && defined(MAX_चळत_ALLOC) && MAX_चळत_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  जर(trans && चळत_alloc_size)
    स_रखो(buffer, 0, MIN(BUFFER_SIZE, आकार(FLOAT) * buffer_size));
#endif

#ifdef SMP

#जर defined(_WIN64) && defined(_M_ARM64)
  जर (m*n > 25000000L) 
    nthreads = अंक_cpu_avail(4);
  अन्यथा
    nthreads = 1;
#अन्यथा
  जर (1L * m * n < 1024L * GEMM_MULTITHREAD_THRESHOLD)
    nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(2);
#endif  

  जर (nthreads == 1) {
#endif  

    (gemv[(पूर्णांक)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अन्यथा {
    (gemv_thread[(पूर्णांक)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, buffer, nthreads);
  }
#endif


  चळत_FREE(buffer);

  FUNCTION_PROFILE_END(4, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  परत;
}
