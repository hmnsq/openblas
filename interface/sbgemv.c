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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#അന്തര്ഭാവ "l1param.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef BGEMM
#define GEMV_THREAD_N bgemv_thread_n
#define GEMV_THREAD_T bgemv_thread_t
#define ERROR_NAME "BGEMV "
#അല്ലെങ്കിലെങ്കിൽ defined(HFLOAT16)
#define GEMV_THREAD_N shgemv_thread_n
#define GEMV_THREAD_T shgemv_thread_t
#define ERROR_NAME "SHGEMV "
#അല്ലെങ്കിൽ
#define GEMV_THREAD_N sbgemv_thread_n
#define GEMV_THREAD_T sbgemv_thread_t
#define ERROR_NAME "SBGEMV "
#endif 

#ifdef SMP
സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, പൂർണ്ണസംഖ്യ) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    അക്ഷരം trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#ifdef SMP
    പൂർണ്ണസംഖ്യ nthreads;
#endif

    പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    എങ്കിൽ (trans == 'N') {i = 0;}
    എങ്കിൽ (trans == 'T') {i = 1;}
    എങ്കിൽ (trans == 'R') {i = 0;}
    എങ്കിൽ (trans == 'C') {i = 1;}

    എങ്കിൽ (incy == 0)       {info = 11;}
    എങ്കിൽ (incx == 0)       {info = 8;}
    എങ്കിൽ (lda < MAX(1, m)) {info = 6;}
    എങ്കിൽ (n < 0)           {info = 3;}
    എങ്കിൽ (m < 0)           {info = 2;}
    എങ്കിൽ (i < 0)           {info = 1;}

    trans = i;

    എങ്കിൽ (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
        തിരികെ;
    }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    പൂർണ്ണസംഖ്യ     trans;
    blasint info,  t;
#ifdef SMP
    പൂർണ്ണസംഖ്യ     nthreads;
#endif

    പൂർണ്ണസംഖ്യ (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    എങ്കിൽ (order == CblasColMajor) {   // Column Major
        എങ്കിൽ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } അല്ലെങ്കിൽ എങ്കിൽ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } അല്ലെങ്കിൽ {                        // Row Major
        എങ്കിൽ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } അല്ലെങ്കിൽ എങ്കിൽ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    എങ്കിൽ (incy == 0)       {info = 11;}
    എങ്കിൽ (incx == 0)       {info = 8;}
    എങ്കിൽ (lda < MAX(1, m)) {info = 6;}
    എങ്കിൽ (n < 0)           {info = 3;}
    എങ്കിൽ (m < 0)           {info = 2;}
    എങ്കിൽ (trans < 0)       {info = 1;}

    എങ്കിൽ (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
        തിരികെ;
    }

#endif

    എങ്കിൽ ((m==0) || (n==0)) തിരികെ;

    എങ്കിൽ (trans) {
        lenx = m;
        leny = n;
    } അല്ലെങ്കിൽ {
        lenx = n;
        leny = m;
    }

#ifdef BGEMM
  ദശാംശം alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#അല്ലെങ്കിൽ
  ദശാംശം alpha_float = alpha;
  ദശാംശം beta_float = beta;
#endif

    എങ്കിൽ (alpha_float == ZERO) {
        എങ്കിൽ (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ശൂന്യ, 0, ശൂന്യ, 0);
        തിരികെ;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    എങ്കിൽ (incx < 0) {x -= (lenx - 1) * incx;}
    എങ്കിൽ (incy < 0) {y -= (leny - 1) * incy;}

#ifdef SMP
    എങ്കിൽ ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    അല്ലെങ്കിൽ
      nthreads = അക്കം_cpu_avail(2);

    എങ്കിൽ (nthreads == 1) {
#endif
        (gemv[(പൂർണ്ണസംഖ്യ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#ifdef SMP
    } അല്ലെങ്കിൽ {
        (gemv_thread[(പൂർണ്ണസംഖ്യ)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#endif

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    തിരികെ;
}
