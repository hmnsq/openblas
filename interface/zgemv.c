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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMV "
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#లేదావేళ defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#లేకపోతే
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  అక్షరం trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *బఫరు;
  పూర్ణాంకం బఫరు_కొలత;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ఒకవేళ (trans == 'N')  i = 0;
  ఒకవేళ (trans == 'T')  i = 1;
  ఒకవేళ (trans == 'R')  i = 2;
  ఒకవేళ (trans == 'C')  i = 3;
  ఒకవేళ (trans == 'O')  i = 4;
  ఒకవేళ (trans == 'U')  i = 5;
  ఒకవేళ (trans == 'S')  i = 6;
  ఒకవేళ (trans == 'D')  i = 7;

  ఒకవేళ (incy == 0)      info = 11;
  ఒకవేళ (incx == 0)      info = 8;
  ఒకవేళ (lda < MAX(1,m)) info = 6;
  ఒకవేళ (n < 0) 	      info = 3;
  ఒకవేళ (m < 0) 	      info = 2;
  ఒకవేళ (i < 0)          info = 1;

  trans = i;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   గణన CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   శూన్య *VALPHA,
	   శూన్య  *va, blasint lda,
	   శూన్య  *vx, blasint incx,
	   శూన్య *VBETA,
	   శూన్య  *vy, blasint incy){

  FLOAT *ALPHA = (FLOAT*) VALPHA;
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  FLOAT *BETA = (FLOAT*) VBETA;
  FLOAT *y = (FLOAT*) vy;
  FLOAT *బఫరు;
  blasint    lenx, leny;
  పూర్ణాంకం trans, బఫరు_కొలత;
  blasint info, t;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (TransA == CblasNoTrans)     trans = 0;
    ఒకవేళ (TransA == CblasTrans)       trans = 1;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 2;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 3;

    info = -1;

    ఒకవేళ (incy == 0)	  info = 11;
    ఒకవేళ (incx == 0)	  info = 8;
    ఒకవేళ (lda < MAX(1, m))  info = 6;
    ఒకవేళ (n < 0)		  info = 3;
    ఒకవేళ (m < 0)		  info = 2;
    ఒకవేళ (trans < 0)        info = 1;

  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (TransA == CblasNoTrans)     trans = 1;
    ఒకవేళ (TransA == CblasTrans)       trans = 0;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 3;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 2;

    info = -1;

    t = n;
    n = m;
    m = t;

    ఒకవేళ (incy == 0)	  info = 11;
    ఒకవేళ (incx == 0)	  info = 8;
    ఒకవేళ (lda < MAX(1, m))  info = 6;
    ఒకవేళ (n < 0)		  info = 3;
    ఒకవేళ (m < 0)		  info = 2;
    ఒకవేళ (trans < 0)        info = 1;

  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#అయితే_ముగింపు

  /*  Quick return if possible. */

  ఒకవేళ (m == 0 || n == 0) తిరిగి;

  lenx = n;
  leny = m;

  ఒకవేళ (trans & 1) lenx = m;
  ఒకవేళ (trans & 1) leny = n;

  ఒకవేళ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ఖాళీ, 0, ఖాళీ, 0);

  ఒకవేళ (alpha_r == ZERO && alpha_i == ZERO) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (lenx - 1) * incx * 2;
  ఒకవేళ (incy < 0) y -= (leny - 1) * incy * 2;

  బఫరు_కొలత = 2 * (m + n) + 128 / పరిమాణం(FLOAT);
#నిర్వచితమైతే WINDOWS_ABI
  బఫరు_కొలత += 160 / పరిమాణం(FLOAT) ;
#అయితే_ముగింపు
  // for alignment
  బఫరు_కొలత = (బఫరు_కొలత + 3) & ~3;
  దొంతర_ALLOC(బఫరు_కొలత, FLOAT, బఫరు);

#ఒకవేళ defined(ARCH_X86_64) && defined(MAX_దొంతర_ALLOC) && MAX_దొంతర_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  ఒకవేళ(trans && దొంతర_alloc_కొలత)
    memset(బఫరు, 0, MIN(బఫరు_కొలత_2, పరిమాణం(FLOAT) * బఫరు_కొలత));
#అయితే_ముగింపు

#నిర్వచితమైతే SMP

#ఒకవేళ defined(_WIN64) && defined(_M_ARM64)
  ఒకవేళ (m*n > 25000000L) 
    nthreads = అంకె_cpu_avail(4);
  లేకపోతే
    nthreads = 1;
#లేకపోతే
  ఒకవేళ (1L * m * n < 1024L * GEMM_MULTITHREAD_పరిమితి)
    nthreads = 1;
  లేకపోతే
    nthreads = అంకె_cpu_avail(2);
#అయితే_ముగింపు  

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు  

    (gemv[(పూర్ణాంకం)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {
    (gemv_thread[(పూర్ణాంకం)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, బఫరు, nthreads);
  }
#అయితే_ముగింపు


  దొంతర_FREE(బఫరు);

  FUNCTION_PROFILE_చివర(4, m * n + m + n,  2 * m * n);

  IDEBUG_చివర;

  తిరిగి;
}
