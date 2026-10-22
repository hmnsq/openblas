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

#समावेश <stdio.h>
#समावेश "common.h"
#समावेश "l1param.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMV "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMV "
#अन्यथा
#घोषणा ERROR_NAME "SGEMV "
#पूर्ण_अगर

#अगर_घोषित SMP
अटल पूर्णांक (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#या_अगर defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#अन्यथा
  sgemv_thread_n, sgemv_thread_t,
#पूर्ण_अगर
};
#पूर्ण_अगर

#अगर_घोषित SMP
#अगर_घोषित DYNAMIC_ARCH
 बाह्य अक्षर* gotoblas_corename(व्योम);
#पूर्ण_अगर

#अगर defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
अटल अंतरभूत पूर्णांक get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, पूर्णांक ncpu) {
  #अगर_घोषित DOUBLE
      वापस (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #अन्यथा
      वापस
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #पूर्ण_अगर
}
#पूर्ण_अगर

#अगर defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
अटल अंतरभूत पूर्णांक get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, पूर्णांक ncpu) {
  वापस
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#पूर्ण_अगर

अटल अंतरभूत पूर्णांक get_gemv_optimal_nthreads(BLASLONG MN) {
  पूर्णांक ncpu = अंक_cpu_avail(3);
#अगर defined(_WIN64) && defined(_M_ARM64)
  अगर (MN > 100000000L)
    वापस अंक_cpu_avail(4);
  वापस 1;
#पूर्ण_अगर
#अगर defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  वापस get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#या_अगर defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  वापस get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#या_अगर defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  अगर (म_भेद(gotoblas_corename(), "neoversev1") == 0) {
    वापस get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  अगर (म_भेद(gotoblas_corename(), "neoversev2") == 0) {
    वापस get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#पूर्ण_अगर

  अगर ( MN < 115200L * GEMM_MULTITHREAD_THRESHOLD )
    वापस 1;
  अन्यथा
    वापस अंक_cpu_avail(2);
}
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *बफर;
  पूर्णांक बफर_आकार;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  अगर (trans == 'N') i = 0;
  अगर (trans == 'T') i = 1;
  अगर (trans == 'R') i = 0;
  अगर (trans == 'C') i = 1;

  अगर (incy == 0)	info = 11;
  अगर (incx == 0)	info = 8;
  अगर (lda < MAX(1, m))	info = 6;
  अगर (n < 0)		info = 3;
  अगर (m < 0)		info = 2;
  अगर (i < 0)          info = 1;

  trans = i;

  अगर (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *बफर;
  blasint lenx, leny;
  पूर्णांक trans, बफर_आकार;
  blasint info, t;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  पूर्णांक (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (TransA == CblasNoTrans)     trans = 0;
    अगर (TransA == CblasTrans)       trans = 1;
    अगर (TransA == CblasConjNoTrans) trans = 0;
    अगर (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    अगर (incy == 0)	  info = 11;
    अगर (incx == 0)	  info = 8;
    अगर (lda < MAX(1, m))  info = 6;
    अगर (n < 0)		  info = 3;
    अगर (m < 0)		  info = 2;
    अगर (trans < 0)        info = 1;

  }

  अगर (order == CblasRowMajor) {
    अगर (TransA == CblasNoTrans)     trans = 1;
    अगर (TransA == CblasTrans)       trans = 0;
    अगर (TransA == CblasConjNoTrans) trans = 1;
    अगर (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    अगर (incy == 0)	  info = 11;
    अगर (incx == 0)	  info = 8;
    अगर (lda < MAX(1, m))  info = 6;
    अगर (n < 0)		  info = 3;
    अगर (m < 0)		  info = 2;
    अगर (trans < 0)        info = 1;

  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर
  अगर ((m==0) || (n==0)) वापस;

  lenx = n;
  leny = m;
  अगर (trans) lenx = m;
  अगर (trans) leny = n;

  अगर (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), शून्य, 0, शून्य, 0);

  अगर (alpha == ZERO) वापस;
	
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (lenx - 1) * incx;
  अगर (incy < 0) y -= (leny - 1) * incy;

  बफर_आकार = m + n + 128 / माप(FLOAT);
#अगर_घोषित WINDOWS_ABI
  बफर_आकार += 160 / माप(FLOAT) ;
#पूर्ण_अगर
  // for alignment
  बफर_आकार = (बफर_आकार + 3) & ~3;
  ढेर_ALLOC(बफर_आकार, FLOAT, बफर);

#अगर_घोषित SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  अगर (nthreads == 1) {
#पूर्ण_अगर

    (gemv[(पूर्णांक)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, बफर);

#अगर_घोषित SMP
  } अन्यथा {

    (gemv_thread[(पूर्णांक)trans])(m, n, alpha, a, lda, x, incx, y, incy, बफर, nthreads);

  }
#पूर्ण_अगर

  ढेर_FREE(बफर);
  FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);

  IDEBUG_END;

  वापस;

}
