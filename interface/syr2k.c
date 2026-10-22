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
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QSYR2K"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DSYR2K"
#अन्यथा
#घोषणा ERROR_NAME "SSYR2K"
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित HEMM
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XSYR2K"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZSYR2K"
#अन्यथा
#घोषणा ERROR_NAME "CSYR2K"
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XHER2K"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHER2K"
#अन्यथा
#घोषणा ERROR_NAME "CHER2K"
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

अटल पूर्णांक (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#अगर_अघोषित HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#अन्यथा
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#पूर्ण_अगर
};

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, अक्षर *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर_घोषित SMP
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

  blasint info;
  पूर्णांक uplo;
  पूर्णांक trans;
  पूर्णांक nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (व्योम *)a;
  args.b = (व्योम *)b;
  args.c = (व्योम *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
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

  अगर (args.ldc < MAX(1,args.n)) info = 12;
  अगर (args.ldb < MAX(1,nrowa))  info =  9;
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
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#अन्यथा
	   व्योम *valpha,
	   व्योम *va, blasint lda,
	   व्योम *vb, blasint ldb,
#पूर्ण_अगर
#अगर !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   व्योम *vbeta,
#पूर्ण_अगर
#अगर_अघोषित COMPLEX
	   FLOAT *c,
#अन्यथा
	   व्योम *vc,
#पूर्ण_अगर
           blasint ldc) {

#अगर_घोषित COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#अगर !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#पूर्ण_अगर
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#पूर्ण_अगर

  blas_arg_t args;
  पूर्णांक uplo, trans;
  blasint info, nrowa;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर_घोषित HEMM
  FLOAT CAlpha[2];
#पूर्ण_अगर

#अगर_घोषित SMP
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

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (व्योम *)a;
  args.b = (व्योम *)b;
  args.c = (व्योम *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#अगर_अघोषित COMPLEX
  args.alpha = (व्योम *)&alpha;
#अन्यथा
  args.alpha = (व्योम *)alpha;
#पूर्ण_अगर

#अगर !defined(COMPLEX) || defined(HEMM)
  args.beta  = (व्योम *)&beta;
#अन्यथा
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

    अगर (args.ldc < MAX(1,args.n)) info = 12;
    अगर (args.ldb < MAX(1,nrowa))  info =  9;
    अगर (args.lda < MAX(1,nrowa))  info =  7;
    अगर (args.k < 0)               info =  4;
    अगर (args.n < 0)               info =  3;
    अगर (trans < 0)                info =  2;
    अगर (uplo  < 0)                info =  1;
  }

  अगर (order == CblasRowMajor) {

#अगर_घोषित HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (व्योम *)CAlpha;
#पूर्ण_अगर

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

    अगर (args.ldc < MAX(1,args.n)) info = 12;
    अगर (args.ldb < MAX(1,nrowa))  info =  9;
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
#अगर defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर
   अगर (order == CblasRowMajor && n == ldc) {
     अगर (Trans == CblasNoTrans && k == lda && k == ldb) {
      अगर (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        वापस;
      }अन्यथा अगर (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        वापस;
      }
     }
     अन्यथा अगर (Trans == CblasTrans && n == lda && n ==ldb) {
      अगर (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        वापस;
      }अन्यथा अगर (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        वापस;
      }
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
  अगर (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = शून्य;
  अगर (args.n*args.k <1000)
	  args.nthreads =1 ;
  अन्यथा
  args.nthreads = अंक_cpu_avail(3);

  अगर (args.nthreads == 1) {
#पूर्ण_अगर

    (syr2k[(uplo << 1) | trans ])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP

  } अन्यथा {

    syrk_thread(mode, &args, शून्य, शून्य, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  वापस;
}
