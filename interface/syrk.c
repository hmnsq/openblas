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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_अपरिभाषितम् COMPLEX
#परिभाषय SMP_देहली_MIN 109944.
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QSYRK "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DSYRK "
#अन्यथा
#परिभाषय ERROR_NAME "SSYRK "
#यदि_अन्तः
#अन्यथा
#परिभाषय SMP_देहली_MIN 14824.
#यदि_अपरिभाषितम् HEMM
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XSYRK "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZSYRK "
#अन्यथा
#परिभाषय ERROR_NAME "CSYRK "
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XHERK "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZHERK "
#अन्यथा
#परिभाषय ERROR_NAME "CHERK "
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि_अपरिभाषितम् GEMM_MULTITHREAD_देहली
#परिभाषय GEMM_MULTITHREAD_देहली 4
#यदि_अन्तः

स्थिर पूर्णाङ्क (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#यदि_अपरिभाषितम् HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#यदि_अन्तः
#अन्यथा
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#यदि_अन्तः
#यदि_अन्तः
};

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;

  blas_arg_t तर्काः;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
  द्विगुण NNK;
#यदि_परिभाषितम् USE_SIMPLE_THREADED_LEVEL3
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क trans;
  पूर्णाङ्क nrowa;

  PRINT_DEBUG_NAME;

  तर्काः.n = *N;
  तर्काः.k = *K;

  तर्काः.a = (शून्य *)a;
  तर्काः.c = (शून्य *)c;

  तर्काः.lda = *ldA;
  तर्काः.ldc = *ldC;

  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;


#यदि_अपरिभाषितम् COMPLEX
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
#यदि_परिभाषितम् HEMM
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
#यदि_अन्तः

#यदि_अन्तः

  nrowa = तर्काः.n;
  यदि (trans & 1) nrowa = तर्काः.k;

  info = 0;

  यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 10;
  यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
  यदि (तर्काः.k < 0)               info =  4;
  यदि (तर्काः.n < 0)               info =  3;
  यदि (trans < 0)                info =  2;
  यदि (uplo  < 0)                info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo, गणना CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#यदि !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#अन्यथा
	   शून्य *valpha,
#यदि_अन्तः
#यदि !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#अन्यथा
	   शून्य *va, blasint lda,
#यदि_अन्तः
#यदि !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   शून्य *vbeta,
#यदि_अन्तः
#यदि !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   शून्य *vc, blasint ldc) {
#यदि_अन्तः

#यदि_परिभाषितम् COMPLEX
#यदि !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#यदि_अन्तः
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#यदि_अन्तः

  blas_arg_t तर्काः;
  पूर्णाङ्क uplo, trans;
  blasint info, nrowa;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
द्विगुण NNK;

#यदि_परिभाषितम् USE_SIMPLE_THREADED_LEVEL3
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  तर्काः.n = n;
  तर्काः.k = k;

  तर्काः.a = (शून्य *)a;
  तर्काः.c = (शून्य *)c;

  तर्काः.lda = lda;
  तर्काः.ldc = ldc;

#यदि !defined(COMPLEX) || defined(HEMM)
  तर्काः.alpha = (शून्य *)&alpha;
  तर्काः.beta  = (शून्य *)&beta;
#अन्यथा
  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;
#यदि_अन्तः

  trans = -1;
  uplo  = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 1;

    यदि (Trans == CblasNoTrans)     trans = 0;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasTrans)       trans = 1;
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 1;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 1;
#यदि_अन्तः

    info = -1;

    nrowa = तर्काः.n;
    यदि (trans & 1) nrowa = तर्काः.k;

    यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 10;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
    यदि (तर्काः.k < 0)               info =  4;
    यदि (तर्काः.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 0;

    यदि (Trans == CblasNoTrans)     trans = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasTrans)       trans = 0;
    यदि (Trans == CblasConjNoTrans) trans = 1;
    यदि (Trans == CblasConjTrans)   trans = 0;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 0;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 0;
#यदि_अन्तः

    info = -1;

    nrowa = तर्काः.n;
    यदि (trans & 1) nrowa = तर्काः.k;

    यदि (तर्काः.ldc < MAX(1,तर्काः.n)) info = 10;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  7;
    यदि (तर्काः.k < 0)               info =  4;
    यदि (तर्काः.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (तर्काः.n == 0) प्रत्यागम;

#यदि_परिभाषितम् DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#यदि_अन्तः
	
#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः
  यदि (क्रमः == CblasRowMajor && n == ldc) {
    यदि (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_ऊर्ध्वस्थम् ? SSYRK_प्रत्यक्षम्_ALPHA_BETA_UN : SSYRK_प्रत्यक्षम्_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      प्रत्यागम;
    } अन्यथा यदि (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_ऊर्ध्वस्थम् ? SSYRK_प्रत्यक्षम्_ALPHA_BETA_UT : SSYRK_प्रत्यक्षम्_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      प्रत्यागम;
    }
  }
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  बफर = (FLOAT *)blas_स्मृतिः_alloc(0);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);

#यदि_परिभाषितम् SMP
#यदि_परिभाषितम् USE_SIMPLE_THREADED_LEVEL3
  यदि (!trans){
    बहुलकम् |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    बहुलकम् |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  बहुलकम् |= (uplo  << BLAS_UPLO_सरणम्);
#यदि_अन्तः

  तर्काः.common = रिक्त;

  NNK = (द्विगुण)(तर्काः.n+1)*(द्विगुण)तर्काः.n*(द्विगुण)तर्काः.k;
  यदि (NNK <= (SMP_देहली_MIN * GEMM_MULTITHREAD_देहली)) {
  तर्काः.nthreads = 1;
  } अन्यथा {
  तर्काः.nthreads = अङ्कः_cpu_avail(3);
  }

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    (syrk[(uplo << 1) | trans ])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_अपरिभाषितम् USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

    syrk_thread(बहुलकम्, &तर्काः, रिक्त, रिक्त, syrk[(uplo << 1) | trans ], sa, sb, तर्काः.nthreads);

#यदि_अन्तः

  }
#यदि_अन्तः

 blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, तर्काः.n * तर्काः.k + तर्काः.n * तर्काः.n / 2, तर्काः.n * तर्काः.n * तर्काः.k);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
