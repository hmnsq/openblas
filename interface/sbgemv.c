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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਅਨ੍ਤਰ੍ਭਾਵ "l1param.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ BGEMM
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_N bgemv_thread_n
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_T bgemv_thread_t
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "BGEMV "
#ਜਾਂਜੇ defined(HFLOAT16)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_N shgemv_thread_n
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_T shgemv_thread_t
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SHGEMV "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_N sbgemv_thread_n
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMV_THREAD_T sbgemv_thread_t
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SBGEMV "
#ਜੇ_ਅੰਤ 

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਸਥਿਰ ਪੂਰਨ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, ਪੂਰਨ) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    ਅੱਖਰ trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
    ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

    ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    ਜੇ (trans == 'N') {i = 0;}
    ਜੇ (trans == 'T') {i = 1;}
    ਜੇ (trans == 'R') {i = 0;}
    ਜੇ (trans == 'C') {i = 1;}

    ਜੇ (incy == 0)       {info = 11;}
    ਜੇ (incx == 0)       {info = 8;}
    ਜੇ (lda < MAX(1, m)) {info = 6;}
    ਜੇ (n < 0)           {info = 3;}
    ਜੇ (m < 0)           {info = 2;}
    ਜੇ (i < 0)           {info = 1;}

    trans = i;

    ਜੇ (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
        ਵਾਪਸ;
    }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    ਪੂਰਨ     trans;
    blasint info,  t;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
    ਪੂਰਨ     nthreads;
#ਜੇ_ਅੰਤ

    ਪੂਰਨ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {   // Column Major
        ਜੇ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } ਨਹੀਂਤਾਂ ਜੇ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } ਨਹੀਂਤਾਂ {                        // Row Major
        ਜੇ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } ਨਹੀਂਤਾਂ ਜੇ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    ਜੇ (incy == 0)       {info = 11;}
    ਜੇ (incx == 0)       {info = 8;}
    ਜੇ (lda < MAX(1, m)) {info = 6;}
    ਜੇ (n < 0)           {info = 3;}
    ਜੇ (m < 0)           {info = 2;}
    ਜੇ (trans < 0)       {info = 1;}

    ਜੇ (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
        ਵਾਪਸ;
    }

#ਜੇ_ਅੰਤ

    ਜੇ ((m==0) || (n==0)) ਵਾਪਸ;

    ਜੇ (trans) {
        lenx = m;
        leny = n;
    } ਨਹੀਂਤਾਂ {
        lenx = n;
        leny = m;
    }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ BGEMM
  ਦਸ਼ਮਲਵ alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#ਨਹੀਂਤਾਂ
  ਦਸ਼ਮਲਵ alpha_float = alpha;
  ਦਸ਼ਮਲਵ beta_float = beta;
#ਜੇ_ਅੰਤ

    ਜੇ (alpha_float == ZERO) {
        ਜੇ (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ਨਿਰਾ, 0, ਨਿਰਾ, 0);
        ਵਾਪਸ;
    }

    IDEBUG_ਸ਼ੁਰੂ;
    FUNCTION_PROFILE_ਸ਼ੁਰੂ();

    ਜੇ (incx < 0) {x -= (lenx - 1) * incx;}
    ਜੇ (incy < 0) {y -= (leny - 1) * incy;}

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
    ਜੇ ( 1L * m * n < 115200L * GEMM_MULTITHREAD_ਹੱਦ )
      nthreads = 1;
    ਨਹੀਂਤਾਂ
      nthreads = ਅੰਕ_cpu_avail(2);

    ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ
        (gemv[(ਪੂਰਨ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
    } ਨਹੀਂਤਾਂ {
        (gemv_thread[(ਪੂਰਨ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#ਜੇ_ਅੰਤ

    FUNCTION_PROFILE_ਸਿਰਾ(1, m * n + m + n,  2 * m * n);
    IDEBUG_ਸਿਰਾ;

    ਵਾਪਸ;
}
