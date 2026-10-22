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

#ضمن <stdio.h>
#ضمن "common.h"
#ضمن "l1param.h"
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_تعریف_شوی BGEMM
#تعریف_کړه GEMV_THREAD_N bgemv_thread_n
#تعریف_کړه GEMV_THREAD_T bgemv_thread_t
#تعریف_کړه ERROR_NAME "BGEMV "
#وإلاإذا defined(HFLOAT16)
#تعریف_کړه GEMV_THREAD_N shgemv_thread_n
#تعریف_کړه GEMV_THREAD_T shgemv_thread_t
#تعریف_کړه ERROR_NAME "SHGEMV "
#وإلا
#تعریف_کړه GEMV_THREAD_N sbgemv_thread_n
#تعریف_کړه GEMV_THREAD_T sbgemv_thread_t
#تعریف_کړه ERROR_NAME "SBGEMV "
#د_که_پای 

#که_تعریف_شوی SMP
ساكن صحيح (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, صحيح) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    حرف trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#که_تعریف_شوی SMP
    صحيح nthreads;
#د_که_پای

    صحيح (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    إذا (trans == 'N') {i = 0;}
    إذا (trans == 'T') {i = 1;}
    إذا (trans == 'R') {i = 0;}
    إذا (trans == 'C') {i = 1;}

    إذا (incy == 0)       {info = 11;}
    إذا (incx == 0)       {info = 8;}
    إذا (lda < MAX(1, m)) {info = 6;}
    إذا (n < 0)           {info = 3;}
    إذا (m < 0)           {info = 2;}
    إذا (i < 0)           {info = 1;}

    trans = i;

    إذا (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
        أرجع;
    }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    صحيح     trans;
    blasint info,  t;
#که_تعریف_شوی SMP
    صحيح     nthreads;
#د_که_پای

    صحيح (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    إذا (ترتیب == Cblas_ستون_Major) {   // Column Major
        إذا (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } وإلا إذا (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } وإلا {                        // Row Major
        إذا (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } وإلا إذا (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    إذا (incy == 0)       {info = 11;}
    إذا (incx == 0)       {info = 8;}
    إذا (lda < MAX(1, m)) {info = 6;}
    إذا (n < 0)           {info = 3;}
    إذا (m < 0)           {info = 2;}
    إذا (trans < 0)       {info = 1;}

    إذا (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
        أرجع;
    }

#د_که_پای

    إذا ((m==0) || (n==0)) أرجع;

    إذا (trans) {
        lenx = m;
        leny = n;
    } وإلا {
        lenx = n;
        leny = m;
    }

#که_تعریف_شوی BGEMM
  عائم alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#وإلا
  عائم alpha_float = alpha;
  عائم beta_float = beta;
#د_که_پای

    إذا (alpha_float == ZERO) {
        إذا (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), عدم, 0, عدم, 0);
        أرجع;
    }

    IDEBUG_پیل;
    FUNCTION_PROFILE_پیل();

    إذا (incx < 0) {x -= (lenx - 1) * incx;}
    إذا (incy < 0) {y -= (leny - 1) * incy;}

#که_تعریف_شوی SMP
    إذا ( 1L * m * n < 115200L * GEMM_MULTITHREAD_پوله )
      nthreads = 1;
    وإلا
      nthreads = عدد_cpu_avail(2);

    إذا (nthreads == 1) {
#د_که_پای
        (gemv[(صحيح)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#که_تعریف_شوی SMP
    } وإلا {
        (gemv_thread[(صحيح)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#د_که_پای

    FUNCTION_PROFILE_پای_ټکی(1, m * n + m + n,  2 * m * n);
    IDEBUG_پای_ټکی;

    أرجع;
}
