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

#include <stdio.h>
#include "common.h"
#include "l1param.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#elif defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#නැතහොත්
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
static පූර්ණ_සංඛ්‍යා (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, පූර්ණ_සංඛ්‍යා) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    var_na trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    පූර්ණ_සංඛ්‍යා nthreads;
#endif

    පූර්ණ_සංඛ්‍යා (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    නම් (trans == 'N') {i = 0;}
    නම් (trans == 'T') {i = 1;}
    නම් (trans == 'R') {i = 0;}
    නම් (trans == 'C') {i = 1;}

    නම් (incy == 0)       {info = 11;}
    නම් (incx == 0)       {info = 8;}
    නම් (lda < MAX(1, m)) {info = 6;}
    නම් (n < 0)           {info = 3;}
    නම් (m < 0)           {info = 2;}
    නම් (i < 0)           {info = 1;}

    trans = i;

    නම් (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
        ආපසු;
    }

#නැතහොත්

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    පූර්ණ_සංඛ්‍යා     trans;
    blasint info,  t;
#ifdef SMP
    පූර්ණ_සංඛ්‍යා     nthreads;
#endif

    පූර්ණ_සංඛ්‍යා (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    නම් (order == CblasColMajor) {   // Column Major
        නම් (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } නැතහොත් නම් (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } නැතහොත් {                        // Row Major
        නම් (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } නැතහොත් නම් (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    නම් (incy == 0)       {info = 11;}
    නම් (incx == 0)       {info = 8;}
    නම් (lda < MAX(1, m)) {info = 6;}
    නම් (n < 0)           {info = 3;}
    නම් (m < 0)           {info = 2;}
    නම් (trans < 0)       {info = 1;}

    නම් (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
        ආපසු;
    }

#endif

    නම් ((m==0) || (n==0)) ආපසු;

    නම් (trans) {
        lenx = m;
        leny = n;
    } නැතහොත් {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  float alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#නැතහොත්
  float alpha_float = alpha;
  float beta_float = beta;
#endif

    නම් (alpha_float == ZERO) {
        නම් (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);
        ආපසු;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    නම් (incx < 0) {x -= (lenx - 1) * incx;}
    නම් (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    නම් ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    නැතහොත්
      nthreads = num_cpu_avail(2);

    නම් (nthreads == 1) {
#endif
        (gemv[(පූර්ණ_සංඛ්‍යා)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } නැතහොත් {
        (gemv_thread[(පූර්ණ_සංඛ්‍යා)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    ආපසු;
}
