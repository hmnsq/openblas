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
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#अन्यथा
#define ERROR_NAME "SSYR2K"
#endif
#अन्यथा
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#अन्यथा
#define ERROR_NAME "CSYR2K"
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#अन्यथा
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

स्थिर पूर्णांक (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#अन्यथा
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *UPLO, अक्षर *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  पूर्णांक uplo;
  पूर्णांक trans;
  पूर्णांक nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (शून्य *)a;
  args.b = (शून्य *)b;
  args.c = (शून्य *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  जर (uplo_arg  == 'U') uplo  = 0;
  जर (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  जर (trans_arg == 'N') trans = 0;
  जर (trans_arg == 'T') trans = 1;
  जर (trans_arg == 'C') trans = 1;
#अन्यथा
#ifdef HEMM
  जर (trans_arg == 'N') trans = 0;
  जर (trans_arg == 'C') trans = 1;
#अन्यथा
  जर (trans_arg == 'N') trans = 0;
  जर (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  जर (trans & 1) nrowa = args.k;

  info = 0;

  जर (args.ldc < MAX(1,args.n)) info = 12;
  जर (args.ldb < MAX(1,nrowa))  info =  9;
  जर (args.lda < MAX(1,nrowa))  info =  7;
  जर (args.k < 0)               info =  4;
  जर (args.n < 0)               info =  3;
  जर (trans < 0)                info =  2;
  जर (uplo  < 0)                info =  1;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, गणना CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
#endif
#जर !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   शून्य *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#अन्यथा
	   शून्य *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#जर !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  पूर्णांक uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (शून्य *)a;
  args.b = (शून्य *)b;
  args.c = (शून्य *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (शून्य *)&alpha;
#अन्यथा
  args.alpha = (शून्य *)alpha;
#endif

#जर !defined(COMPLEX) || defined(HEMM)
  args.beta  = (शून्य *)&beta;
#अन्यथा
  args.beta  = (शून्य *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  जर (order == CblasColMajor) {
    जर (Uplo == CblasUpper) uplo  = 0;
    जर (Uplo == CblasLower) uplo  = 1;

    जर (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    जर (Trans == CblasTrans)       trans = 1;
    जर (Trans == CblasConjNoTrans) trans = 0;
    जर (Trans == CblasConjTrans)   trans = 1;
#किंवाजर !defined(HEMM)
    जर (Trans == CblasTrans)       trans = 1;
#अन्यथा
    जर (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    जर (trans & 1) nrowa = args.k;

    जर (args.ldc < MAX(1,args.n)) info = 12;
    जर (args.ldb < MAX(1,nrowa))  info =  9;
    जर (args.lda < MAX(1,nrowa))  info =  7;
    जर (args.k < 0)               info =  4;
    जर (args.n < 0)               info =  3;
    जर (trans < 0)                info =  2;
    जर (uplo  < 0)                info =  1;
  }

  जर (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (शून्य *)CAlpha;
#endif

    जर (Uplo == CblasUpper) uplo  = 1;
    जर (Uplo == CblasLower) uplo  = 0;

    जर (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    जर (Trans == CblasTrans)       trans = 0;
    जर (Trans == CblasConjNoTrans) trans = 1;
    जर (Trans == CblasConjTrans)   trans = 0;
#किंवाजर !defined(HEMM)
    जर (Trans == CblasTrans)       trans = 0;
#अन्यथा
    जर (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    जर (trans & 1) nrowa = args.k;

    जर (args.ldc < MAX(1,args.n)) info = 12;
    जर (args.ldb < MAX(1,nrowa))  info =  9;
    जर (args.lda < MAX(1,nrowa))  info =  7;
    जर (args.k < 0)               info =  4;
    जर (args.n < 0)               info =  3;
    जर (trans < 0)                info =  2;
    जर (uplo  < 0)                info =  1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

  जर (args.n == 0) परत;

#ifdef DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(शून्य);
#endif

#जर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#जर defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
जर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#जर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   जर (order == CblasRowMajor && n == ldc) {
     जर (Trans == CblasNoTrans && k == lda && k == ldb) {
      जर (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        परत;
      }अन्यथा जर (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        परत;
      }
     }
     अन्यथा जर (Trans == CblasTrans && n == lda && n ==ldb) {
      जर (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        परत;
      }अन्यथा जर (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        परत;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  जर (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = रिक्त;
  जर (args.n*args.k <1000)
	  args.nthreads =1 ;
  अन्यथा
  args.nthreads = अंक_cpu_avail(3);

  जर (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP

  } अन्यथा {

    syrk_thread(mode, &args, रिक्त, रिक्त, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  परत;
}
