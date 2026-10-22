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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#અન્તર્ભાવ "l1param.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#અથવાજો defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#નહીંતર
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
સ્થિર પૂર્ણાંક (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, પૂર્ણાંક) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    અક્ષર trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    પૂર્ણાંક nthreads;
#endif

    પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    જો (trans == 'N') {i = 0;}
    જો (trans == 'T') {i = 1;}
    જો (trans == 'R') {i = 0;}
    જો (trans == 'C') {i = 1;}

    જો (incy == 0)       {info = 11;}
    જો (incx == 0)       {info = 8;}
    જો (lda < MAX(1, m)) {info = 6;}
    જો (n < 0)           {info = 3;}
    જો (m < 0)           {info = 2;}
    જો (i < 0)           {info = 1;}

    trans = i;

    જો (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
        પાછા;
    }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    પૂર્ણાંક     trans;
    blasint info,  t;
#ifdef SMP
    પૂર્ણાંક     nthreads;
#endif

    પૂર્ણાંક (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    જો (order == CblasColMajor) {   // Column Major
        જો (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } નહીંતર જો (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } નહીંતર {                        // Row Major
        જો (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } નહીંતર જો (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    જો (incy == 0)       {info = 11;}
    જો (incx == 0)       {info = 8;}
    જો (lda < MAX(1, m)) {info = 6;}
    જો (n < 0)           {info = 3;}
    જો (m < 0)           {info = 2;}
    જો (trans < 0)       {info = 1;}

    જો (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
        પાછા;
    }

#endif

    જો ((m==0) || (n==0)) પાછા;

    જો (trans) {
        lenx = m;
        leny = n;
    } નહીંતર {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  દશાંશ alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#નહીંતર
  દશાંશ alpha_float = alpha;
  દશાંશ beta_float = beta;
#endif

    જો (alpha_float == ZERO) {
        જો (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ખાલી, 0, ખાલી, 0);
        પાછા;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    જો (incx < 0) {x -= (lenx - 1) * incx;}
    જો (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    જો ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    નહીંતર
      nthreads = અંક_cpu_avail(2);

    જો (nthreads == 1) {
#endif
        (gemv[(પૂર્ણાંક)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } નહીંતર {
        (gemv_thread[(પૂર્ણાંક)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    પાછા;
}
