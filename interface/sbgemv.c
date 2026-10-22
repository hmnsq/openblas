/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project.                              */
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

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ BGEMM
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_N bgemv_thread_n
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_T bgemv_thread_t
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "BGEMV "
#ಅಥವಾಆದರೆ defined(HFLOAT16)
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_N shgemv_thread_n
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_T shgemv_thread_t
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SHGEMV "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_N sbgemv_thread_n
#ವ್ಯಾಖ್ಯಾನಿಸು GEMV_THREAD_T sbgemv_thread_t
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SBGEMV "
#ಆದರೆ_ಅಂತ್ಯ 

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, ಪೂರ್ಣಾಂಕ) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    ಅಕ್ಷರ trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
    ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

    ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    ಆದರೆ (trans == 'N') {i = 0;}
    ಆದರೆ (trans == 'T') {i = 1;}
    ಆದರೆ (trans == 'R') {i = 0;}
    ಆದರೆ (trans == 'C') {i = 1;}

    ಆದರೆ (incy == 0)       {info = 11;}
    ಆದರೆ (incx == 0)       {info = 8;}
    ಆದರೆ (lda < MAX(1, m)) {info = 6;}
    ಆದರೆ (n < 0)           {info = 3;}
    ಆದರೆ (m < 0)           {info = 2;}
    ಆದರೆ (i < 0)           {info = 1;}

    trans = i;

    ಆದರೆ (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
        ಹಿಂತಿರುಗಿ;
    }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    ಪೂರ್ಣಾಂಕ     trans;
    blasint info,  t;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
    ಪೂರ್ಣಾಂಕ     nthreads;
#ಆದರೆ_ಅಂತ್ಯ

    ಪೂರ್ಣಾಂಕ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {   // Column Major
        ಆದರೆ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } ಇಲ್ಲದಿದ್ದರೆ {                        // Row Major
        ಆದರೆ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    ಆದರೆ (incy == 0)       {info = 11;}
    ಆದರೆ (incx == 0)       {info = 8;}
    ಆದರೆ (lda < MAX(1, m)) {info = 6;}
    ಆದರೆ (n < 0)           {info = 3;}
    ಆದರೆ (m < 0)           {info = 2;}
    ಆದರೆ (trans < 0)       {info = 1;}

    ಆದರೆ (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
        ಹಿಂತಿರುಗಿ;
    }

#ಆದರೆ_ಅಂತ್ಯ

    ಆದರೆ ((m==0) || (n==0)) ಹಿಂತಿರುಗಿ;

    ಆದರೆ (trans) {
        lenx = m;
        leny = n;
    } ಇಲ್ಲದಿದ್ದರೆ {
        lenx = n;
        leny = m;
    }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ BGEMM
  ದಶಮಾಂಶ alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#ಇಲ್ಲದಿದ್ದರೆ
  ದಶಮಾಂಶ alpha_float = alpha;
  ದಶಮಾಂಶ beta_float = beta;
#ಆದರೆ_ಅಂತ್ಯ

    ಆದರೆ (alpha_float == ZERO) {
        ಆದರೆ (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ಖಾಲಿ, 0, ಖಾಲಿ, 0);
        ಹಿಂತಿರುಗಿ;
    }

    IDEBUG_ಆರಂಭ;
    FUNCTION_PROFILE_ಆರಂಭ();

    ಆದರೆ (incx < 0) {x -= (lenx - 1) * incx;}
    ಆದರೆ (incy < 0) {y -= (leny - 1) * incy;}

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
    ಆದರೆ ( 1L * m * n < 115200L * GEMM_MULTITHREAD_ಮಿತಿ )
      nthreads = 1;
    ಇಲ್ಲದಿದ್ದರೆ
      nthreads = ಅಂಕೆ_cpu_avail(2);

    ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ
        (gemv[(ಪೂರ್ಣಾಂಕ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
    } ಇಲ್ಲದಿದ್ದರೆ {
        (gemv_thread[(ಪೂರ್ಣಾಂಕ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#ಆದರೆ_ಅಂತ್ಯ

    FUNCTION_PROFILE_ಕೊನೆ(1, m * n + m + n,  2 * m * n);
    IDEBUG_ಕೊನೆ;

    ಹಿಂತಿರುಗಿ;
}
