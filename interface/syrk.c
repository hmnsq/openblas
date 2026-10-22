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

#समावेश <stdio.h>
#समावेश <ctype.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_अघोषित COMPLEX
#घोषणा SMP_THRESHOLD_MIN 109944.
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QSYRK "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DSYRK "
#अन्यथा
#घोषणा ERROR_NAME "SSYRK "
#पूर्ण_अगर
#अन्यथा
#घोषणा SMP_THRESHOLD_MIN 14824.
#अगर_अघोषित HEMM
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XSYRK "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZSYRK "
#अन्यथा
#घोषणा ERROR_NAME "CSYRK "
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XHERK "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHERK "
#अन्यथा
#घोषणा ERROR_NAME "CHERK "
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर_अघोषित GEMM_MULTITHREAD_THRESHOLD
#घोषणा GEMM_MULTITHREAD_THRESHOLD 4
#पूर्ण_अगर

अटल पूर्णांक (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#अगर_अघोषित HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#पूर्ण_अगर
#अन्यथा
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#पूर्ण_अगर
#पूर्ण_अगर
};

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, अक्षर *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर_घोषित SMP
  द्विगुन NNK;
#अगर_घोषित USE_SIMPLE_THREADED_LEVEL3
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

  blasint info;
  पूर्णांक uplo;
  पूर्णांक trans;
  पूर्णांक nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (व्योम *)a;
  args.c = (व्योम *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;


#अगर_अघोषित COMPLEX
  अगर (trans_arg == 'N') trans = 0;
  अगर (trans_arg == 'T') trans = 1;
  अगर (trans_arg == 'C') trans = 1;
#अन्यथा
#अगर_घोषित HEMM
  अगर (trans_arg == 'N') trans = 0;
  अगर (trans_arg == 'C') trans = 1;
#अन्यथा
  अगर (trans_arg == 'N') trans = 0;
  अगर (trans_arg == 'T') trans = 1;
#पूर्ण_अगर

#पूर्ण_अगर

  nrowa = args.n;
  अगर (trans & 1) nrowa = args.k;

  info = 0;

  अगर (args.ldc < MAX(1,args.n)) info = 10;
  अगर (args.lda < MAX(1,nrowa))  info =  7;
  अगर (args.k < 0)               info =  4;
  अगर (args.n < 0)               info =  3;
  अगर (trans < 0)                info =  2;
  अगर (uplo  < 0)                info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_UPLO Uplo, क्रमागत CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#अगर !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#अन्यथा
	   व्योम *valpha,
#पूर्ण_अगर
#अगर !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#अन्यथा
	   व्योम *va, blasint lda,
#पूर्ण_अगर
#अगर !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   व्योम *vbeta,
#पूर्ण_अगर
#अगर !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   व्योम *vc, blasint ldc) {
#पूर्ण_अगर

#अगर_घोषित COMPLEX
#अगर !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#पूर्ण_अगर
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#पूर्ण_अगर

  blas_arg_t args;
  पूर्णांक uplo, trans;
  blasint info, nrowa;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर_घोषित SMP
द्विगुन NNK;

#अगर_घोषित USE_SIMPLE_THREADED_LEVEL3
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (व्योम *)a;
  args.c = (व्योम *)c;

  args.lda = lda;
  args.ldc = ldc;

#अगर !defined(COMPLEX) || defined(HEMM)
  args.alpha = (व्योम *)&alpha;
  args.beta  = (व्योम *)&beta;
#अन्यथा
  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;
#पूर्ण_अगर

  trans = -1;
  uplo  = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (Uplo == CblasUpper) uplo  = 0;
    अगर (Uplo == CblasLower) uplo  = 1;

    अगर (Trans == CblasNoTrans)     trans = 0;
#अगर_अघोषित COMPLEX
    अगर (Trans == CblasTrans)       trans = 1;
    अगर (Trans == CblasConjNoTrans) trans = 0;
    अगर (Trans == CblasConjTrans)   trans = 1;
#या_अगर !defined(HEMM)
    अगर (Trans == CblasTrans)       trans = 1;
#अन्यथा
    अगर (Trans == CblasConjTrans)   trans = 1;
#पूर्ण_अगर

    info = -1;

    nrowa = args.n;
    अगर (trans & 1) nrowa = args.k;

    अगर (args.ldc < MAX(1,args.n)) info = 10;
    अगर (args.lda < MAX(1,nrowa))  info =  7;
    अगर (args.k < 0)               info =  4;
    अगर (args.n < 0)               info =  3;
    अगर (trans < 0)                info =  2;
    अगर (uplo  < 0)                info =  1;
  }

  अगर (order == CblasRowMajor) {
    अगर (Uplo == CblasUpper) uplo  = 1;
    अगर (Uplo == CblasLower) uplo  = 0;

    अगर (Trans == CblasNoTrans)     trans = 1;
#अगर_अघोषित COMPLEX
    अगर (Trans == CblasTrans)       trans = 0;
    अगर (Trans == CblasConjNoTrans) trans = 1;
    अगर (Trans == CblasConjTrans)   trans = 0;
#या_अगर !defined(HEMM)
    अगर (Trans == CblasTrans)       trans = 0;
#अन्यथा
    अगर (Trans == CblasConjTrans)   trans = 0;
#पूर्ण_अगर

    info = -1;

    nrowa = args.n;
    अगर (trans & 1) nrowa = args.k;

    अगर (args.ldc < MAX(1,args.n)) info = 10;
    अगर (args.lda < MAX(1,nrowa))  info =  7;
    अगर (args.k < 0)               info =  4;
    अगर (args.n < 0)               info =  3;
    अगर (trans < 0)                info =  2;
    अगर (uplo  < 0)                info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

  अगर (args.n == 0) वापस;

#अगर_घोषित DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(व्योम);
#पूर्ण_अगर
	
#अगर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#अगर defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर
  अगर (order == CblasRowMajor && n == ldc) {
    अगर (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      वापस;
    } अन्यथा अगर (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      वापस;
    }
  }
#पूर्ण_अगर
#पूर्ण_अगर

#पूर्ण_अगर


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  बफर = (FLOAT *)blas_memory_alloc(0);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#अगर_घोषित SMP
#अगर_घोषित USE_SIMPLE_THREADED_LEVEL3
  अगर (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#पूर्ण_अगर

  args.common = शून्य;

  NNK = (द्विगुन)(args.n+1)*(द्विगुन)args.n*(द्विगुन)args.k;
  अगर (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } अन्यथा {
  args.nthreads = अंक_cpu_avail(3);
  }

  अगर (args.nthreads == 1) {
#पूर्ण_अगर

    (syrk[(uplo << 1) | trans ])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_अघोषित USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, शून्य, शून्य, sa, sb, 0);

#अन्यथा

    syrk_thread(mode, &args, शून्य, शून्य, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#पूर्ण_अगर

  }
#पूर्ण_अगर

 blas_memory_free(बफर);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  वापस;
}
