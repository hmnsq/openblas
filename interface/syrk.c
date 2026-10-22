/*********************************************************************/
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
#అన్తర్భావ <ctype.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు SMP_పరిమితి_MIN 109944.
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QSYRK "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DSYRK "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SSYRK "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచించు SMP_పరిమితి_MIN 14824.
#నిర్వచితం_కాకపోతే HEMM
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XSYRK "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZSYRK "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CSYRK "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XHERK "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZHERK "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CHERK "
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే GEMM_MULTITHREAD_పరిమితి
#నిర్వచించు GEMM_MULTITHREAD_పరిమితి 4
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#నిర్వచితం_కాకపోతే HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#అయితే_ముగింపు
#లేకపోతే
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#అయితే_ముగింపు
#అయితే_ముగింపు
};

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *UPLO, అక్షరం *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  అక్షరం uplo_arg  = *UPLO;
  అక్షరం trans_arg = *TRANS;

  blas_arg_t వాదనలు;

  FLOAT *బఫరు;
  FLOAT *sa, *sb;

#నిర్వచితమైతే SMP
  ద్విగుణ NNK;
#నిర్వచితమైతే USE_SIMPLE_THREADED_LEVEL3
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

  blasint info;
  పూర్ణాంకం uplo;
  పూర్ణాంకం trans;
  పూర్ణాంకం nrowa;

  PRINT_DEBUG_NAME;

  వాదనలు.n = *N;
  వాదనలు.k = *K;

  వాదనలు.a = (శూన్య *)a;
  వాదనలు.c = (శూన్య *)c;

  వాదనలు.lda = *ldA;
  వాదనలు.ldc = *ldC;

  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;


#నిర్వచితం_కాకపోతే COMPLEX
  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
  ఒకవేళ (trans_arg == 'C') trans = 1;
#లేకపోతే
#నిర్వచితమైతే HEMM
  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'C') trans = 1;
#లేకపోతే
  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
#అయితే_ముగింపు

#అయితే_ముగింపు

  nrowa = వాదనలు.n;
  ఒకవేళ (trans & 1) nrowa = వాదనలు.k;

  info = 0;

  ఒకవేళ (వాదనలు.ldc < MAX(1,వాదనలు.n)) info = 10;
  ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  7;
  ఒకవేళ (వాదనలు.k < 0)               info =  4;
  ఒకవేళ (వాదనలు.n < 0)               info =  3;
  ఒకవేళ (trans < 0)                info =  2;
  ఒకవేళ (uplo  < 0)                info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_UPLO Uplo, గణన CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ఒకవేళ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#లేకపోతే
	   శూన్య *valpha,
#అయితే_ముగింపు
#ఒకవేళ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#లేకపోతే
	   శూన్య *va, blasint lda,
#అయితే_ముగింపు
#ఒకవేళ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#లేకపోతే
	   శూన్య *vbeta,
#అయితే_ముగింపు
#ఒకవేళ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#లేకపోతే
	   శూన్య *vc, blasint ldc) {
#అయితే_ముగింపు

#నిర్వచితమైతే COMPLEX
#ఒకవేళ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#అయితే_ముగింపు
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#అయితే_ముగింపు

  blas_arg_t వాదనలు;
  పూర్ణాంకం uplo, trans;
  blasint info, nrowa;

  FLOAT *బఫరు;
  FLOAT *sa, *sb;

#నిర్వచితమైతే SMP
ద్విగుణ NNK;

#నిర్వచితమైతే USE_SIMPLE_THREADED_LEVEL3
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  వాదనలు.n = n;
  వాదనలు.k = k;

  వాదనలు.a = (శూన్య *)a;
  వాదనలు.c = (శూన్య *)c;

  వాదనలు.lda = lda;
  వాదనలు.ldc = ldc;

#ఒకవేళ !defined(COMPLEX) || defined(HEMM)
  వాదనలు.alpha = (శూన్య *)&alpha;
  వాదనలు.beta  = (శూన్య *)&beta;
#లేకపోతే
  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;
#అయితే_ముగింపు

  trans = -1;
  uplo  = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (Uplo == Cblas_పై) uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 1;

    ఒకవేళ (Trans == CblasNoTrans)     trans = 0;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (Trans == CblasTrans)       trans = 1;
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 0;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 1;
#లేదావేళ !defined(HEMM)
    ఒకవేళ (Trans == CblasTrans)       trans = 1;
#లేకపోతే
    ఒకవేళ (Trans == CblasConjTrans)   trans = 1;
#అయితే_ముగింపు

    info = -1;

    nrowa = వాదనలు.n;
    ఒకవేళ (trans & 1) nrowa = వాదనలు.k;

    ఒకవేళ (వాదనలు.ldc < MAX(1,వాదనలు.n)) info = 10;
    ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  7;
    ఒకవేళ (వాదనలు.k < 0)               info =  4;
    ఒకవేళ (వాదనలు.n < 0)               info =  3;
    ఒకవేళ (trans < 0)                info =  2;
    ఒకవేళ (uplo  < 0)                info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (Uplo == Cblas_పై) uplo  = 1;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 0;

    ఒకవేళ (Trans == CblasNoTrans)     trans = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (Trans == CblasTrans)       trans = 0;
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 1;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 0;
#లేదావేళ !defined(HEMM)
    ఒకవేళ (Trans == CblasTrans)       trans = 0;
#లేకపోతే
    ఒకవేళ (Trans == CblasConjTrans)   trans = 0;
#అయితే_ముగింపు

    info = -1;

    nrowa = వాదనలు.n;
    ఒకవేళ (trans & 1) nrowa = వాదనలు.k;

    ఒకవేళ (వాదనలు.ldc < MAX(1,వాదనలు.n)) info = 10;
    ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  7;
    ఒకవేళ (వాదనలు.k < 0)               info =  4;
    ఒకవేళ (వాదనలు.n < 0)               info =  3;
    ఒకవేళ (trans < 0)                info =  2;
    ఒకవేళ (uplo  < 0)                info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

  ఒకవేళ (వాదనలు.n == 0) తిరిగి;

#నిర్వచితమైతే DYNAMIC_ARCH
బాహ్య అక్షరం* gotoblas_corename(శూన్య);
#అయితే_ముగింపు
	
#ఒకవేళ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ఒకవేళ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
ఒకవేళ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ఒకవేళ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#అయితే_ముగింపు
)
#అయితే_ముగింపు
  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && n == ldc) {
    ఒకవేళ (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_పై ? SSYRK_ప్రత్యక్ష_ALPHA_BETA_UN : SSYRK_ప్రత్యక్ష_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      తిరిగి;
    } లేకపోతే ఒకవేళ (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_పై ? SSYRK_ప్రత్యక్ష_ALPHA_BETA_UT : SSYRK_ప్రత్యక్ష_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      తిరిగి;
    }
  }
#అయితే_ముగింపు
#అయితే_ముగింపు

#అయితే_ముగింపు


  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(0);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

  sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);

#నిర్వచితమైతే SMP
#నిర్వచితమైతే USE_SIMPLE_THREADED_LEVEL3
  ఒకవేళ (!trans){
    బహుళకం |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } లేకపోతే {
    బహుళకం |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  బహుళకం |= (uplo  << BLAS_UPLO_జరుపు);
#అయితే_ముగింపు

  వాదనలు.common = ఖాళీ;

  NNK = (ద్విగుణ)(వాదనలు.n+1)*(ద్విగుణ)వాదనలు.n*(ద్విగుణ)వాదనలు.k;
  ఒకవేళ (NNK <= (SMP_పరిమితి_MIN * GEMM_MULTITHREAD_పరిమితి)) {
  వాదనలు.nthreads = 1;
  } లేకపోతే {
  వాదనలు.nthreads = అంకె_cpu_avail(3);
  }

  ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    (syrk[(uplo << 1) | trans ])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP
  } లేకపోతే {

#నిర్వచితం_కాకపోతే USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#లేకపోతే

    syrk_thread(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, syrk[(uplo << 1) | trans ], sa, sb, వాదనలు.nthreads);

#అయితే_ముగింపు

  }
#అయితే_ముగింపు

 blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE, వాదనలు.n * వాదనలు.k + వాదనలు.n * వాదనలు.n / 2, వాదనలు.n * వాదనలు.n * వాదనలు.k);

  IDEBUG_చివర;

  తిరిగి;
}
