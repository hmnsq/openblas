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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#అన్తర్భావ "l1param.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితమైతే BGEMM
#నిర్వచించు GEMV_THREAD_N bgemv_thread_n
#నిర్వచించు GEMV_THREAD_T bgemv_thread_t
#నిర్వచించు ERROR_NAME "BGEMV "
#లేదావేళ defined(HFLOAT16)
#నిర్వచించు GEMV_THREAD_N shgemv_thread_n
#నిర్వచించు GEMV_THREAD_T shgemv_thread_t
#నిర్వచించు ERROR_NAME "SHGEMV "
#లేకపోతే
#నిర్వచించు GEMV_THREAD_N sbgemv_thread_n
#నిర్వచించు GEMV_THREAD_T sbgemv_thread_t
#నిర్వచించు ERROR_NAME "SBGEMV "
#అయితే_ముగింపు 

#నిర్వచితమైతే SMP
స్థిర పూర్ణాంకం (*gemv_thread[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG, పూర్ణాంకం) = {
    GEMV_THREAD_N, GEMV_THREAD_T,
};
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *TRANS, blasint *M, blasint *N, FLOAT *ALPHA, IFLOAT *a, blasint *LDA, IFLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{
    అక్షరం trans = *TRANS;
    blasint m = *M;
    blasint n = *N;
    blasint lda = *LDA;
    blasint incx = *INCX;
    blasint incy = *INCY;
    FLOAT alpha = *ALPHA;
    FLOAT beta  = *BETA;
#నిర్వచితమైతే SMP
    పూర్ణాంకం nthreads;
#అయితే_ముగింపు

    పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG, IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    blasint info;
    blasint lenx, leny;
    blasint i;

    PRINT_DEBUG_NAME;

    TOUPPER(trans);

    info = 0;

    i = -1;

    ఒకవేళ (trans == 'N') {i = 0;}
    ఒకవేళ (trans == 'T') {i = 1;}
    ఒకవేళ (trans == 'R') {i = 0;}
    ఒకవేళ (trans == 'C') {i = 1;}

    ఒకవేళ (incy == 0)       {info = 11;}
    ఒకవేళ (incx == 0)       {info = 8;}
    ఒకవేళ (lda < MAX(1, m)) {info = 6;}
    ఒకవేళ (n < 0)           {info = 3;}
    ఒకవేళ (m < 0)           {info = 2;}
    ఒకవేళ (i < 0)           {info = 1;}

    trans = i;

    ఒకవేళ (info != 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
        తిరిగి;
    }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_TRANSPOSE TransA, blasint m, blasint n, FLOAT alpha, IFLOAT *a, blasint lda, IFLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{
    blasint lenx,  leny;
    పూర్ణాంకం     trans;
    blasint info,  t;
#నిర్వచితమైతే SMP
    పూర్ణాంకం     nthreads;
#అయితే_ముగింపు

    పూర్ణాంకం (*gemv[])(BLASLONG, BLASLONG, FLOAT, IFLOAT *, BLASLONG,  IFLOAT * , BLASLONG, FLOAT, FLOAT *, BLASLONG) = {
        GEMV_N, GEMV_T,
    };

    PRINT_DEBUG_CNAME;

    trans = -1;
    info  =  0;

    ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {   // Column Major
        ఒకవేళ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 0;
        } లేకపోతే ఒకవేళ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 1;
        }
    } లేకపోతే {                        // Row Major
        ఒకవేళ (TransA == CblasNoTrans || TransA == CblasConjNoTrans) {
            trans = 1;
        } లేకపోతే ఒకవేళ (TransA == CblasTrans || TransA == CblasConjTrans) {
            trans = 0;
        }

        t = n;
        n = m;
        m = t;
    }

    info = -1;

    ఒకవేళ (incy == 0)       {info = 11;}
    ఒకవేళ (incx == 0)       {info = 8;}
    ఒకవేళ (lda < MAX(1, m)) {info = 6;}
    ఒకవేళ (n < 0)           {info = 3;}
    ఒకవేళ (m < 0)           {info = 2;}
    ఒకవేళ (trans < 0)       {info = 1;}

    ఒకవేళ (info >= 0) {
        BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
        తిరిగి;
    }

#అయితే_ముగింపు

    ఒకవేళ ((m==0) || (n==0)) తిరిగి;

    ఒకవేళ (trans) {
        lenx = m;
        leny = n;
    } లేకపోతే {
        lenx = n;
        leny = m;
    }

#నిర్వచితమైతే BGEMM
  దశాంశం alpha_float, beta_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
  SBF16TOS_K(1, &beta, 1, &beta_float, 1);
#లేకపోతే
  దశాంశం alpha_float = alpha;
  దశాంశం beta_float = beta;
#అయితే_ముగింపు

    ఒకవేళ (alpha_float == ZERO) {
        ఒకవేళ (beta_float != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), ఖాళీ, 0, ఖాళీ, 0);
        తిరిగి;
    }

    IDEBUG_ప్రారంభం;
    FUNCTION_PROFILE_ప్రారంభం();

    ఒకవేళ (incx < 0) {x -= (lenx - 1) * incx;}
    ఒకవేళ (incy < 0) {y -= (leny - 1) * incy;}

#నిర్వచితమైతే SMP
    ఒకవేళ ( 1L * m * n < 115200L * GEMM_MULTITHREAD_పరిమితి )
      nthreads = 1;
    లేకపోతే
      nthreads = అంకె_cpu_avail(2);

    ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు
        (gemv[(పూర్ణాంకం)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy);
#నిర్వచితమైతే SMP
    } లేకపోతే {
        (gemv_thread[(పూర్ణాంకం)trans])(m, n, alpha, a, lda, x, incx, beta, y, incy, nthreads);
    }
#అయితే_ముగింపు

    FUNCTION_PROFILE_చివర(1, m * n + m + n,  2 * m * n);
    IDEBUG_చివర;

    తిరిగి;
}
