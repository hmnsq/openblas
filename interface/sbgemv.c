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

#समावेश <stdio.h>
#समावेश "common.h"
#समावेश "l1param.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित BGEMM
#घोषणा GEMV_THREAD_N bgemv_thread_n
#घोषणा GEMV_THREAD_T bgemv_thread_t
#घोषणा ERROR_NAME "BGEMV "
#या_अगर defined(HFLOAT16)
#घोषणा GEMV_THREAD_N shgemv_thread_n
#घोषणा GEMV_THREAD_T shgemv_thread_t
#घोषणा ERROR_NAME "SHGEMV "
#अन्यथा
#घोषणा GEMV_THREAD_N sbgemv_thread_n
#घोषणा GEMV_THREAD_T sbgemv_thread_t
#घोषणा ERROR_NAME "SBGEMV "
#पूर्ण_अगर 

#अगर_घोषित SMP
अटल पूर्णांक (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, पूर्णांक) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    अक्षर trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#अगर_घोषित SMP
    पूर्णांक nthreads;
#पूर्ण_अगर

    पूर्णांक (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    अगर (trans == 'N') {i = 0;}
    अगर (trans == 'T') {i = 1;}
    अगर (trans == 'R') {i = 0;}
    अगर (trans == 'C') {i = 1;}

    अगर (incy == 0)       {info = 11;}
    अगर (incx == 0)       {info = 8;}
    अगर (lda < MAX(1, m)) {info = 6;}
    अगर (n < 0)           {info = 3;}
    अगर (m < 0)           {info = 2;}
    अगर (i < 0)           {info = 1;}

    trans = i;

    अगर (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
        वापस;
    }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    पूर्णांक     trans;
    blasint info,  t;
#अगर_घोषित SMP
    पूर्णांक     nthreads;
#पूर्ण_अगर

    पूर्णांक (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    अगर (order == CblasColMajor) {   // Column Major
        अगर (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } अन्यथा अगर (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } अन्यथा {                        // Row Major
        अगर (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } अन्यथा अगर (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    अगर (incy == 0)       {info = 11;}
    अगर (incx == 0)       {info = 8;}
    अगर (lda < MAX(1, m)) {info = 6;}
    अगर (n < 0)           {info = 3;}
    अगर (m < 0)           {info = 2;}
    अगर (trans < 0)       {info = 1;}

    अगर (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
        वापस;
    }

#पूर्ण_अगर

    अगर ((m==0) || (n==0)) वापस;

    अगर (trans) {
        lenx = m;
        leny = n;
    } अन्यथा {
        lenx = n;
        leny = m;
    }

#अगर_घोषित BGEMM
  भग्न alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#अन्यथा
  भग्न alpha_float = alpha;
  भग्न beta_float = beta;
#पूर्ण_अगर

    अगर (alpha_float == ZERO) {
        अगर (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), शून्य, 0, शून्य, 0);
        वापस;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    अगर (incx < 0) {x -= (lenx - 1) * incx;}
    अगर (incy < 0) {y -= (leny - 1) * incy;}

#अगर_घोषित SMP
    अगर ( 1L * m * n < 115200L * GEMM_MULTITHREAD_THRESHOLD )
      nthreads = 1;
    अन्यथा
      nthreads = अंक_cpu_avail(2);

    अगर (nthreads == 1) {
#पूर्ण_अगर
        (gemv[(पूर्णांक)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#अगर_घोषित SMP
    } अन्यथा {
        (gemv_thread[(पूर्णांक)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#पूर्ण_अगर

    FUNCTION_PROFILE_END(1, m * n + m + n,  2 * m * n);
    IDEBUG_END;

    वापस;
}
