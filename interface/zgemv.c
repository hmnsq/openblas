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
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMV "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMV "
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT *, FLOAT *, BLASLONG,  FLOAT * , BLASLONG, FLOAT *, BLASLONG, FLOAT *, ಪೂರ್ಣಾಂಕ) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  xgemv_thread_n, xgemv_thread_t, xgemv_thread_r, xgemv_thread_c, xgemv_thread_o, xgemv_thread_u, xgemv_thread_s, xgemv_thread_d,
#ಅಥವಾಆದರೆ defined DOUBLE
  zgemv_thread_n, zgemv_thread_t, zgemv_thread_r, zgemv_thread_c, zgemv_thread_o, zgemv_thread_u, zgemv_thread_s, zgemv_thread_d,
#ಇಲ್ಲದಿದ್ದರೆ
  cgemv_thread_n, cgemv_thread_t, cgemv_thread_r, cgemv_thread_c, cgemv_thread_o, cgemv_thread_u, cgemv_thread_s, cgemv_thread_d,
#ಆದರೆ_ಅಂತ್ಯ
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANS, blasint *M, blasint *N,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA,  FLOAT *y, blasint *INCY){

  ಅಕ್ಷರ trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;

  FLOAT *ಬಫರ್;
  ಪೂರ್ಣಾಂಕ ಬಫರ್_ಗಾತ್ರ;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ಆದರೆ (trans == 'N')  i = 0;
  ಆದರೆ (trans == 'T')  i = 1;
  ಆದರೆ (trans == 'R')  i = 2;
  ಆದರೆ (trans == 'C')  i = 3;
  ಆದರೆ (trans == 'O')  i = 4;
  ಆದರೆ (trans == 'U')  i = 5;
  ಆದರೆ (trans == 'S')  i = 6;
  ಆದರೆ (trans == 'D')  i = 7;

  ಆದರೆ (incy == 0)      info = 11;
  ಆದರೆ (incx == 0)      info = 8;
  ಆದರೆ (lda < MAX(1,m)) info = 6;
  ಆದರೆ (n < 0) 	      info = 3;
  ಆದರೆ (m < 0) 	      info = 2;
  ಆದರೆ (i < 0)          info = 1;

  trans = i;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE TransA,
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
  FLOAT *ಬಫರ್;
  blasint    lenx, leny;
  ಪೂರ್ಣಾಂಕ trans, ಬಫರ್_ಗಾತ್ರ;
  blasint info, t;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

  ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG,
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

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ಆದರೆ (TransA == CblasNoTrans)     trans = 0;
    ಆದರೆ (TransA == CblasTrans)       trans = 1;
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 2;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 3;

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
    ಆದರೆ (TransA == CblasConjNoTrans) trans = 3;
    ಆದರೆ (TransA == CblasConjTrans)   trans = 2;

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

  /*  Quick return if possible. */

  ಆದರೆ (m == 0 || n == 0) ಹಿಂತಿರುಗಿ;

  lenx = n;
  leny = m;

  ಆದರೆ (trans & 1) lenx = m;
  ಆದರೆ (trans & 1) leny = n;

  ಆದರೆ (beta_r != ONE || beta_i != ZERO) SCAL_K(leny, 0, 0, beta_r, beta_i, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);

  ಆದರೆ (alpha_r == ZERO && alpha_i == ZERO) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (lenx - 1) * incx * 2;
  ಆದರೆ (incy < 0) y -= (leny - 1) * incy * 2;

  ಬಫರ್_ಗಾತ್ರ = 2 * (m + n) + 128 / ಗಾತ್ರ(FLOAT);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ WINDOWS_ABI
  ಬಫರ್_ಗಾತ್ರ += 160 / ಗಾತ್ರ(FLOAT) ;
#ಆದರೆ_ಅಂತ್ಯ
  // for alignment
  ಬಫರ್_ಗಾತ್ರ = (ಬಫರ್_ಗಾತ್ರ + 3) & ~3;
  ರಾಶಿ_ALLOC(ಬಫರ್_ಗಾತ್ರ, FLOAT, ಬಫರ್);

#ಆದರೆ defined(ARCH_X86_64) && defined(MAX_ರಾಶಿ_ALLOC) && MAX_ರಾಶಿ_ALLOC > 0
  // cgemv_t.S return NaN if there are NaN or Inf in the buffer (see bug #746)
  ಆದರೆ(trans && ರಾಶಿ_alloc_ಗಾತ್ರ)
    memset(ಬಫರ್, 0, MIN(ಬಫರ್_ಗಾತ್ರ_2, ಗಾತ್ರ(FLOAT) * ಬಫರ್_ಗಾತ್ರ));
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP

#ಆದರೆ defined(_WIN64) && defined(_M_ARM64)
  ಆದರೆ (m*n > 25000000L) 
    nthreads = ಅಂಕೆ_cpu_avail(4);
  ಇಲ್ಲದಿದ್ದರೆ
    nthreads = 1;
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (1L * m * n < 1024L * GEMM_MULTITHREAD_ಮಿತಿ)
    nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ
    nthreads = ಅಂಕೆ_cpu_avail(2);
#ಆದರೆ_ಅಂತ್ಯ  

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ  

    (gemv[(ಪೂರ್ಣಾಂಕ)trans])(m, n, 0, alpha_r, alpha_i, a, lda, x, incx, y, incy, ಬಫರ್);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {
    (gemv_thread[(ಪೂರ್ಣಾಂಕ)trans])(m, n, ALPHA, a, lda, x, incx, y, incy, ಬಫರ್, nthreads);
  }
#ಆದರೆ_ಅಂತ್ಯ


  ರಾಶಿ_FREE(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(4, m * n + m + n,  2 * m * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
