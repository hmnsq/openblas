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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#অন্তর্ভাৱ "l1param.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#নোচেদ্যদি defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#অন্যথা
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
স্থির পূর্ণাঙ্ক (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, পূর্ণাঙ্ক) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(ৱর্ণ *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    ৱর্ণ trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    পূর্ণাঙ্ক nthreads;
#endif

    পূর্ণাঙ্ক (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    যদি (trans == 'N') {i = 0;}
    যদি (trans == 'T') {i = 1;}
    যদি (trans == 'R') {i = 0;}
    যদি (trans == 'C') {i = 1;}

    যদি (incy == 0)       {info = 11;}
    যদি (incx == 0)       {info = 8;}
    যদি (lda < MAX(1, m)) {info = 6;}
    যদি (n < 0)           {info = 3;}
    যদি (m < 0)           {info = 2;}
    যদি (i < 0)           {info = 1;}

    trans = i;

    যদি (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
        ঘূৰাই;
    }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    পূর্ণাঙ্ক     trans;
    blasint info,  t;
#ifdef SMP
    পূর্ণাঙ্ক     nthreads;
#endif

    পূর্ণাঙ্ক (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    যদি (order == CblasColMajor) {   // Column Major
        যদি (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } অন্যথা যদি (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } অন্যথা {                        // Row Major
        যদি (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } অন্যথা যদি (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    যদি (incy == 0)       {info = 11;}
    যদি (incx == 0)       {info = 8;}
    যদি (lda < MAX(1, m)) {info = 6;}
    যদি (n < 0)           {info = 3;}
    যদি (m < 0)           {info = 2;}
    যদি (trans < 0)       {info = 1;}

    যদি (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
        ঘূৰাই;
    }

#endif

    যদি ((m==0) || (n==0)) ঘূৰাই;

    যদি (trans) {
        lenx = m;
        leny = n;
    } অন্যথা {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  ভিন্ন alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#অন্যথা
  ভিন্ন alpha_float = alpha;
  ভিন্ন beta_float = beta;
#endif

    যদি (alpha_float == ZERO) {
        যদি (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), রিক্ত, 0, রিক্ত, 0);
        ঘূৰাই;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    যদি (incx < 0) {x -= (lenx - 1) * incx;}
    যদি (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    যদি ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    অন্যথা
      nthreads = অঙ্ক_cpu_avail(2);

    যদি (nthreads == 1) {
#endif
        (gemv[(পূর্ণাঙ্ক)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } অন্যথা {
        (gemv_thread[(পূর্ণাঙ্ক)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    ঘূৰাই;
}
