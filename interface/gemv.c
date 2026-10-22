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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ಅನ್ತರ್ಭಾವ "l1param.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGEMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGEMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGEMV "
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  qgemv_thread_n, qgemv_thread_t,
#ಅಥವಾಆದರೆ defined DOUBLE
  dgemv_thread_n, dgemv_thread_t,
#ಇಲ್ಲದಿದ್ದರೆ
  sgemv_thread_n, sgemv_thread_t,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DYNAMIC_ARCH
 ಬಾಹ್ಯ ಅಕ್ಷರ* gotoblas_corename(ಶೂನ್ಯ);
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(DYNAMIC_ARCH) || defined(NEOVERSEV1)
ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemv_optimal_nthreads_neoversev1(BLASLONG MN, ಪೂರ್ಣಾಂಕ ncpu) {
  #ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DOUBLE
      ಹಿಂತಿರುಗಿ (MN < 8100L)      ? 1
          : (MN < 12100L)     ? MIN(ncpu, 2)
          : (MN < 36100L)     ? MIN(ncpu, 4)
          : (MN < 84100L)     ? MIN(ncpu, 8)
          : (MN < 348100L)    ? MIN(ncpu, 16)
          : (MN < 435600L)    ? MIN(ncpu, 24)
          : (MN < 810000L)    ? MIN(ncpu, 32)
          : (MN < 1050625L)   ? MIN(ncpu, 40)
          : ncpu;
  #ಇಲ್ಲದಿದ್ದರೆ
      ಹಿಂತಿರುಗಿ
            (MN < 25600L)     ? 1
          : (MN < 63001L)     ? MIN(ncpu, 4)
          : (MN < 202500L)    ? MIN(ncpu, 8)
          : (MN < 806404L)    ? MIN(ncpu, 16)
          : (MN < 1638400L)    ? MIN(ncpu, 32)
          : ncpu;
  #ಆದರೆ_ಅಂತ್ಯ
}
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(DYNAMIC_ARCH) || defined(NEOVERSEV2)
ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemv_optimal_nthreads_neoversev2(BLASLONG MN, ಪೂರ್ಣಾಂಕ ncpu) {
  ಹಿಂತಿರುಗಿ
      MN < 24964L    ? 1
    : MN < 145924L    ? MIN(ncpu, 8)
    : MN < 692224L   ? MIN(ncpu, 16)
    : MN < 1638400L  ? MIN(ncpu, 32)
    : ncpu;
}
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ ಸಾಲೊಳಗೆ ಪೂರ್ಣಾಂಕ get_gemv_optimal_nthreads(BLASLONG MN) {
  ಪೂರ್ಣಾಂಕ ncpu = ಅಂಕೆ_cpu_avail(3);
#ಆದರೆ defined(_WIN64) && defined(_M_ARM64)
  ಆದರೆ (MN > 100000000L)
    ಹಿಂತಿರುಗಿ ಅಂಕೆ_cpu_avail(4);
  ಹಿಂತಿರುಗಿ 1;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ defined(NEOVERSEV1) && !defined(COMPLEX) && !defined(BFLOAT16)
  ಹಿಂತಿರುಗಿ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
#ಅಥವಾಆದರೆ defined(NEOVERSEV2) && !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)
  ಹಿಂತಿರುಗಿ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
#ಅಥವಾಆದರೆ defined(DYNAMIC_ARCH) && !defined(COMPLEX) && !defined(BFLOAT16)
  ಆದರೆ (strcmp(gotoblas_corename(), "neoversev1") == 0) {
    ಹಿಂತಿರುಗಿ get_gemv_optimal_nthreads_neoversev1(MN, ncpu);
  }
  ಆದರೆ (strcmp(gotoblas_corename(), "neoversev2") == 0) {
    ಹಿಂತಿರುಗಿ get_gemv_optimal_nthreads_neoversev2(MN, ncpu);
  }
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ ( MN < 115200L * GEMM_MULTITHREAD_ಮಿತಿ )
    ಹಿಂತಿರುಗಿ 1;
  ಇಲ್ಲದಿದ್ದರೆ
    ಹಿಂತಿರುಗಿ ಅಂಕೆ_cpu_avail(2);
}
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANS, blasint *M, blasint *N,
	   FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	   FLOAT *x, blasint *INCX,
	   FLOAT *BETA, FLOAT *y, blasint *INCY){

  ಅಕ್ಷರ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;
  FLOAT *ಬಫರ್;
  ಪೂರ್ಣಾಂಕ ಬಫರ್_ಗಾತ್ರ;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  ಆದರೆ (trans == 'N') i = 0;
  ಆದರೆ (trans == 'T') i = 1;
  ಆದರೆ (trans == 'R') i = 0;
  ಆದರೆ (trans == 'C') i = 1;

  ಆದರೆ (incy == 0)	info = 11;
  ಆದರೆ (incx == 0)	info = 8;
  ಆದರೆ (lda < MAX(1, m))	info = 6;
  ಆದರೆ (n < 0)		info = 3;
  ಆದರೆ (m < 0)		info = 2;
  ಆದರೆ (i < 0)          info = 1;

  trans = i;

  ಆದರೆ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *ಬಫರ್;
  blasint lenx, leny;
  ಪೂರ್ಣಾಂಕ trans, ಬಫರ್_ಗಾತ್ರ;
  blasint info, t;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
    GEMV_N, GEMV_T,
  };

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (TransA == CblasNoTrans)     trans = 0;
    ಆದರೆ (TransA == CblasTrans)       trans = 1;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 0;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    ಆದರೆ (incy == 0)	  info = 11;
    ಆದರೆ (incx == 0)	  info = 8;
    ಆದರೆ (lda < MAX(1, m))  info = 6;
    ಆದರೆ (n < 0)		  info = 3;
    ಆದರೆ (m < 0)		  info = 2;
    ಆದರೆ (trans < 0)        info = 1;

  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ಆದರೆ (TransA == CblasNoTrans)     trans = 1;
    ಆದರೆ (TransA == CblasTrans)       trans = 0;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 1;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    ಆದರೆ (incy == 0)	  info = 11;
    ಆದರೆ (incx == 0)	  info = 8;
    ಆದರೆ (lda < MAX(1, m))  info = 6;
    ಆದರೆ (n < 0)		  info = 3;
    ಆದರೆ (m < 0)		  info = 2;
    ಆದರೆ (trans < 0)        info = 1;

  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಆದರೆ_ಅಂತ್ಯ
  ಆದರೆ ((m==0) || (n==0)) ಹಿಂತಿರುಗಿ;

  lenx = n;
  leny = m;
  ಆದರೆ (trans) lenx = m;
  ಆದರೆ (trans) leny = n;

  ಆದರೆ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);

  ಆದರೆ (alpha == ZERO) ಹಿಂತಿರುಗಿ;
	
  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (lenx - 1) * incx;
  ಆದರೆ (incy < 0) y -= (leny - 1) * incy;

  ಬಫರ್_ಗಾತ್ರ = m + n + 128 / ಗಾತ್ರ(FLOAT);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ WINDOWS_ABI
  ಬಫರ್_ಗಾತ್ರ += 160 / ಗಾತ್ರ(FLOAT) ;
#ಆದರೆ_ಅಂತ್ಯ
  // for alignment
  ಬಫರ್_ಗಾತ್ರ = (ಬಫರ್_ಗಾತ್ರ + 3) & ~3;
  ರಾಶಿ_ALLOC(ಬಫರ್_ಗಾತ್ರ, FLOAT, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  nthreads = get_gemv_optimal_nthreads(1L * m * n);

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (gemv[(ಪೂರ್ಣಾಂಕ)trans])(m, n, 0, alpha, a, lda, x, incx, y, incy, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

    (gemv_thread[(ಪೂರ್ಣಾಂಕ)trans])(m, n, alpha, a, lda, x, incx, y, incy, ಬಫರ್, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  ರಾಶಿ_FREE(ಬಫರ್);
  FUNCTION_PROFILE_ಕೊನೆ(1, m * n + m + n,  2 * m * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;

}
