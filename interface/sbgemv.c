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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#अन्तर्भाव "l1param.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् BGEMM
#परिभाषय GEMV_THREAD_N bgemv_thread_n
#परिभाषय GEMV_THREAD_T bgemv_thread_t
#परिभाषय ERROR_NAME "BGEMV "
#नोचेद्यदि defined(HFLOAT16)
#परिभाषय GEMV_THREAD_N shgemv_thread_n
#परिभाषय GEMV_THREAD_T shgemv_thread_t
#परिभाषय ERROR_NAME "SHGEMV "
#अन्यथा
#परिभाषय GEMV_THREAD_N sbgemv_thread_n
#परिभाषय GEMV_THREAD_T sbgemv_thread_t
#परिभाषय ERROR_NAME "SBGEMV "
#यदि_अन्तः 

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, पूर्णाङ्क) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    वर्ण trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#यदि_परिभाषितम् SMP
    पूर्णाङ्क nthreads;
#यदि_अन्तः

    पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    यदि (trans == 'N') {i = 0;}
    यदि (trans == 'T') {i = 1;}
    यदि (trans == 'R') {i = 0;}
    यदि (trans == 'C') {i = 1;}

    यदि (incy == 0)       {info = 11;}
    यदि (incx == 0)       {info = 8;}
    यदि (lda < MAX(1, m)) {info = 6;}
    यदि (n < 0)           {info = 3;}
    यदि (m < 0)           {info = 2;}
    यदि (i < 0)           {info = 1;}

    trans = i;

    यदि (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
        प्रत्यागम;
    }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    पूर्णाङ्क     trans;
    blasint info,  t;
#यदि_परिभाषितम् SMP
    पूर्णाङ्क     nthreads;
#यदि_अन्तः

    पूर्णाङ्क (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    यदि (क्रमः == CblasColMajor) {   // Column Major
        यदि (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } अन्यथा यदि (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } अन्यथा {                        // Row Major
        यदि (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } अन्यथा यदि (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    यदि (incy == 0)       {info = 11;}
    यदि (incx == 0)       {info = 8;}
    यदि (lda < MAX(1, m)) {info = 6;}
    यदि (n < 0)           {info = 3;}
    यदि (m < 0)           {info = 2;}
    यदि (trans < 0)       {info = 1;}

    यदि (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
        प्रत्यागम;
    }

#यदि_अन्तः

    यदि ((m==0) || (n==0)) प्रत्यागम;

    यदि (trans) {
        lenx = m;
        leny = n;
    } अन्यथा {
        lenx = n;
        leny = m;
    }

#यदि_परिभाषितम् BGEMM
  भिन्न alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#अन्यथा
  भिन्न alpha_float = alpha;
  भिन्न beta_float = beta;
#यदि_अन्तः

    यदि (alpha_float == ZERO) {
        यदि (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);
        प्रत्यागम;
    }

    IDEBUG_START;
    FUNCTION_PROFILE_START();

    यदि (incx < 0) {x -= (lenx - 1) * incx;}
    यदि (incy < 0) {y -= (leny - 1) * incy;}

#यदि_परिभाषितम् SMP
    यदि ( 1L * m * n < 115200L * GEMM_MULTITHREAD_देहली )
      nthreads = 1;
    अन्यथा
      nthreads = अङ्कः_cpu_avail(2);

    यदि (nthreads == 1) {
#यदि_अन्तः
        (gemv[(पूर्णाङ्क)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#यदि_परिभाषितम् SMP
    } अन्यथा {
        (gemv_thread[(पूर्णाङ्क)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#यदि_अन्तः

    FUNCTION_PROFILE_अन्तः(1, m * n + m + n,  2 * m * n);
    IDEBUG_अन्तः;

    प्रत्यागम;
}
