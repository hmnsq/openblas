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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#அந்தர்பாவ "l1param.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#அல்லதெனில் defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#இல்லையெனில்
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
நிலையான முழு (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, முழு) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    எழுத்து trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    முழு nthreads;
#endif

    முழு (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    எனில் (trans == 'N') {i = 0;}
    எனில் (trans == 'T') {i = 1;}
    எனில் (trans == 'R') {i = 0;}
    எனில் (trans == 'C') {i = 1;}

    எனில் (incy == 0)       {info = 11;}
    எனில் (incx == 0)       {info = 8;}
    எனில் (lda < MAX(1, m)) {info = 6;}
    எனில் (n < 0)           {info = 3;}
    எனில் (m < 0)           {info = 2;}
    எனில் (i < 0)           {info = 1;}

    trans = i;

    எனில் (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
        திரும்பு;
    }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    முழு     trans;
    blasint info,  t;
#ifdef SMP
    முழு     nthreads;
#endif

    முழு (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    எனில் (order == CblasColMajor) {   // Column Major
        எனில் (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } இல்லையெனில் எனில் (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } இல்லையெனில் {                        // Row Major
        எனில் (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } இல்லையெனில் எனில் (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    எனில் (incy == 0)       {info = 11;}
    எனில் (incx == 0)       {info = 8;}
    எனில் (lda < MAX(1, m)) {info = 6;}
    எனில் (n < 0)           {info = 3;}
    எனில் (m < 0)           {info = 2;}
    எனில் (trans < 0)       {info = 1;}

    எனில் (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
        திரும்பு;
    }

#endif

    எனில் ((m==0) || (n==0)) திரும்பு;

    எனில் (trans) {
        lenx = m;
        leny = n;
    } இல்லையெனில் {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  தசம alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#இல்லையெனில்
  தசம alpha_float = alpha;
  தசம beta_float = beta;
#endif

    எனில் (alpha_float == ZERO) {
        எனில் (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), வெற்று, 0, வெற்று, 0);
        திரும்பு;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    எனில் (incx < 0) {x -= (lenx - 1) * incx;}
    எனில் (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    எனில் ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    இல்லையெனில்
      nthreads = இலக்கம்_cpu_avail(2);

    எனில் (nthreads == 1) {
#endif
        (gemv[(முழு)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } இல்லையெனில் {
        (gemv_thread[(முழு)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    திரும்பு;
}
