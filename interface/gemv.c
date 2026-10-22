/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project                               */
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
#అన్తర్భావ "l1param.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGEMV "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEMV "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEMV "
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, పూర్ణాంకం) = {
#నిర్వచితమైతే XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#లేదావేళ defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#లేకపోతే
  sgemv_thread_n, sgemv_thread_t,
#అయితే_ముగింపు
};
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
#నిర్వచితమైతే DYNAMIC_ARCH
 బాహ్య అక్షరం* gotoblas_corename(శూన్య);
#అయితే_ముగింపు

#ఒకవేళ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
స్థిర అంతఃవరుస పూర్ణాంకం get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, పూర్ణాంకం ncpu) {
  #నిర్వచితమైతే DOUBLE
      తిరిగి (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #లేకపోతే
      తిరిగి
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #అయితే_ముగింపు
}
#అయితే_ముగింపు

#ఒకవేళ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
స్థిర అంతఃవరుస పూర్ణాంకం get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, పూర్ణాంకం ncpu) {
  తిరిగి
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#అయితే_ముగింపు

స్థిర అంతఃవరుస పూర్ణాంకం get_gemv_optimal_nthreads(BLASLONG MN) {
  పూర్ణాంకం ncpu = అంకె_cpu_avail(3);
#ఒకవేళ defined(_WIN64) && defined(_M_ARM64)
  ఒకవేళ (MN > 100000000L)
    తిరిగి అంకె_cpu_avail(4);
  తిరిగి 1;
#అయితే_ముగింపు
#ఒకవేళ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  తిరిగి get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#లేదావేళ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  తిరిగి get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#లేదావేళ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ఒకవేళ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    తిరిగి get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ఒకవేళ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    తిరిగి get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#అయితే_ముగింపు

  ఒకవేళ ( MN < 115200L * GEMM_MULTITHREAD_పరిమితి )
    తిరిగి 1;
  లేకపోతే
    తిరిగి అంకె_cpu_avail(2);
}
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  అక్షరం trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *బఫరు;
  పూర్ణాంకం బఫరు_కొలత;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ఒకవేళ (trans == 'N') i = 0;
  ఒకవేళ (trans == 'T') i = 1;
  ఒకవేళ (trans == 'R') i = 0;
  ఒకవేళ (trans == 'C') i = 1;

  ఒకవేళ (incy == 0)	info = 11;
  ఒకవేళ (incx == 0)	info = 8;
  ఒకవేళ (lda < MAX(1, m))	info = 6;
  ఒకవేళ (n < 0)		info = 3;
  ఒకవేళ (m < 0)		info = 2;
  ఒకవేళ (i < 0)          info = 1;

  trans = i;

  ఒకవేళ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   గణన CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *బఫరు;
  blasint lenx, leny;
  పూర్ణాంకం trans, బఫరు_కొలత;
  blasint info, t;
#నిర్వచితమైతే SMP
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

  పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (TransA == CblasNoTrans)     trans = 0;
    ఒకవేళ (TransA == CblasTrans)       trans = 1;
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 0;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 1;

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
    ఒకవేళ (TransA == CblasConjNoTrans) trans = 1;
    ఒకవేళ (TransA == CblasConjTrans)   trans = 0;

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
  ఒకవేళ ((m==0) || (n==0)) తిరిగి;

  lenx = n;
  leny = m;
  ఒకవేళ (trans) lenx = m;
  ఒకవేళ (trans) leny = n;

  ఒకవేళ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ఖాళీ, 0, ఖాళీ, 0);

  ఒకవేళ (alpha == ZERO) తిరిగి;
	
  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (lenx - 1) * incx;
  ఒకవేళ (incy < 0) y -= (leny - 1) * incy;

  బఫరు_కొలత = m + n + 128 / పరిమాణం(FLOAT);
#నిర్వచితమైతే WINDOWS_ABI
  బఫరు_కొలత += 160 / పరిమాణం(FLOAT) ;
#అయితే_ముగింపు
  // for alignment
  బఫరు_కొలత = (బఫరు_కొలత + 3) & ~3;
  దొంతర_ALLOC(బఫరు_కొలత, FLOAT, బఫరు);

#నిర్వచితమైతే SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

    (gemv[(పూర్ణాంకం)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, బఫరు);

#నిర్వచితమైతే SMP
  } లేకపోతే {

    (gemv_thread[(పూర్ణాంకం)trans])(m, n, alpha, a, lda, x, incx, y, incy, బఫరు, nthreads);

  }
#అయితే_ముగింపు

  దొంతర_FREE(బఫరు);
  FUNCTION_PROFILE_చివర(1, m * n + m + n,  2 * m * n);

  IDEBUG_చివర;

  తిరిగి;

}
