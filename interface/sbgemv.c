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
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#ನೋಚೇದ್ಯದಿ defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#ಅನ್ಯಥಾ
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, ಪೂರ್ಣಾಙ್ಕ) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    ವರ್ಣ trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    ಪೂರ್ಣಾಙ್ಕ nthreads;
#endif

    ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    ಯದಿ (trans == 'N') {i = 0;}
    ಯದಿ (trans == 'T') {i = 1;}
    ಯದಿ (trans == 'R') {i = 0;}
    ಯದಿ (trans == 'C') {i = 1;}

    ಯದಿ (incy == 0)       {info = 11;}
    ಯದಿ (incx == 0)       {info = 8;}
    ಯದಿ (lda < MAX(1, m)) {info = 6;}
    ಯದಿ (n < 0)           {info = 3;}
    ಯದಿ (m < 0)           {info = 2;}
    ಯದಿ (i < 0)           {info = 1;}

    trans = i;

    ಯದಿ (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
        ಪ್ರತ್ಯಾಗಮ;
    }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    ಪೂರ್ಣಾಙ್ಕ     trans;
    blasint info,  t;
#ifdef SMP
    ಪೂರ್ಣಾಙ್ಕ     nthreads;
#endif

    ಪೂರ್ಣಾಙ್ಕ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    ಯದಿ (order == CblasColMajor) {   // Column Major
        ಯದಿ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } ಅನ್ಯಥಾ ಯದಿ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } ಅನ್ಯಥಾ {                        // Row Major
        ಯದಿ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } ಅನ್ಯಥಾ ಯದಿ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    ಯದಿ (incy == 0)       {info = 11;}
    ಯದಿ (incx == 0)       {info = 8;}
    ಯದಿ (lda < MAX(1, m)) {info = 6;}
    ಯದಿ (n < 0)           {info = 3;}
    ಯದಿ (m < 0)           {info = 2;}
    ಯದಿ (trans < 0)       {info = 1;}

    ಯದಿ (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
        ಪ್ರತ್ಯಾಗಮ;
    }

#endif

    ಯದಿ ((m==0) || (n==0)) ಪ್ರತ್ಯಾಗಮ;

    ಯದಿ (trans) {
        lenx = m;
        leny = n;
    } ಅನ್ಯಥಾ {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  ಭಿನ್ನ alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#ಅನ್ಯಥಾ
  ಭಿನ್ನ alpha_float = alpha;
  ಭಿನ್ನ beta_float = beta;
#endif

    ಯದಿ (alpha_float == ZERO) {
        ಯದಿ (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ರಿಕ್ತ, 0, ರಿಕ್ತ, 0);
        ಪ್ರತ್ಯಾಗಮ;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    ಯದಿ (incx < 0) {x -= (lenx - 1) * incx;}
    ಯದಿ (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    ಯದಿ ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    ಅನ್ಯಥಾ
      nthreads = ಅಂಕೆ_cpu_avail(2);

    ಯದಿ (nthreads == 1) {
#endif
        (gemv[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } ಅನ್ಯಥಾ {
        (gemv_thread[(ಪೂರ್ಣಾಙ್ಕ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    ಪ್ರತ್ಯಾಗಮ;
}
