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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMV "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMV "
#अन्यथा
#परिभाषय ERROR_NAME "CGEMV "
#यदि_अन्तः

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#यदि_परिभाषितम् XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#नोचेद्यदि defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#अन्यथा
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#यदि_अन्तः
};
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

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

  FLOAT *बफर;
  पूर्णाङ्क बफर_आकार;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

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

  यदि (trans == 'N')  i = 0;
  यदि (trans == 'T')  i = 1;
  यदि (trans == 'R')  i = 2;
  यदि (trans == 'C')  i = 3;
  यदि (trans == 'O')  i = 4;
  यदि (trans == 'U')  i = 5;
  यदि (trans == 'S')  i = 6;
  यदि (trans == 'D')  i = 7;

  यदि (incy == 0)      info = 11;
  यदि (incx == 0)      info = 8;
  यदि (lda < MAX(1,m)) info = 6;
  यदि (n < 0) 	      info = 3;
  यदि (m < 0) 	      info = 2;
  यदि (i < 0)          info = 1;

  trans = i;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः,
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
  FLOAT *बफर;
  blasint    lenx, leny;
  पूर्णाङ्क trans, बफर_आकार;
  blasint info, t;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

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

  यदि (क्रमः == CblasColMajor) {
    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 2;
    यदि (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    यदि (incy == 0)	  info = 11;
    यदि (incx == 0)	  info = 8;
    यदि (lda < MAX(1, m))  info = 6;
    यदि (n < 0)		  info = 3;
    यदि (m < 0)		  info = 2;
    यदि (trans < 0)        info = 1;

  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 3;
    यदि (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    यदि (incy == 0)	  info = 11;
    यदि (incx == 0)	  info = 8;
    यदि (lda < MAX(1, m))  info = 6;
    यदि (n < 0)		  info = 3;
    यदि (m < 0)		  info = 2;
    यदि (trans < 0)        info = 1;

  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  /*  Quick return if possible. */

  यदि (m == 0 || n == 0) प्रत्यागम;

  lenx = n;
  leny = m;

  यदि (trans & 1) lenx = m;
  यदि (trans & 1) leny = n;

  यदि (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  यदि (alpha_r == ZERO && alpha_i == ZERO) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (lenx - 1) * incx * 2;
  यदि (incy < 0) y -= (leny - 1) * incy * 2;

  बफर_आकार = 2 * (m + n) + 128 / परिमाण(FLOAT);
#यदि_परिभाषितम् WINDOWS_ABI
  बफर_आकार += 160 / परिमाण(FLOAT) ;
#यदि_अन्तः
  // for alignment
  बफर_आकार = (बफर_आकार + 3) & ~3;
  राशिः_ALLOC(बफर_आकार, FLOAT, बफर);

#यदि defined(ARCH_X86_64) && defined(MAX_राशिः_ALLOC) && MAX_राशिः_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  यदि(trans && राशिः_alloc_आकार)
    स्मृतिं_पूरय(बफर, 0, MIN(बफर_आकार_2, परिमाण(FLOAT) * बफर_आकार));
#यदि_अन्तः

#यदि_परिभाषितम् SMP

#यदि defined(_WIN64) && defined(_M_ARM64)
  यदि (m*n > 25000000L) 
    nthreads = अङ्कः_cpu_avail(4);
  अन्यथा
    nthreads = 1;
#अन्यथा
  यदि (1L * m * n < 1024L * GEMM_MULTITHREAD_देहली)
    nthreads = 1;
  अन्यथा
    nthreads = अङ्कः_cpu_avail(2);
#यदि_अन्तः  

  यदि (nthreads == 1) {
#यदि_अन्तः  

    (gemv[(पूर्णाङ्क)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, बफर);

#यदि_परिभाषितम् SMP
  } अन्यथा {
    (gemv_thread[(पूर्णाङ्क)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, बफर, nthreads);
  }
#यदि_अन्तः


  राशिः_FREE(बफर);

  FUNCTION_PROFILE_अन्तः(4, m * n + m + n,  2 * m * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
