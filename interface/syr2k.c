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
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#अण्णहा
#define ERROR_NAME "SSYR2K"
#endif
#अण्णहा
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#अण्णहा
#define ERROR_NAME "CSYR2K"
#endif
#अण्णहा
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#अण्णहा
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

स्थिर पूर्णाङ्क (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#अण्णहा
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अण्णहा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अण्णहा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अण्णहा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क trans;
  पूर्णाङ्क nrowa;

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

  जइ (uplo_arg  == 'U') uplo  = 0;
  जइ (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  जइ (trans_arg == 'N') trans = 0;
  जइ (trans_arg == 'T') trans = 1;
  जइ (trans_arg == 'C') trans = 1;
#अण्णहा
#ifdef HEMM
  जइ (trans_arg == 'N') trans = 0;
  जइ (trans_arg == 'C') trans = 1;
#अण्णहा
  जइ (trans_arg == 'N') trans = 0;
  जइ (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  जइ (trans & 1) nrowa = args.k;

  info = 0;

  जइ (args.ldc < MAX(1,args.n)) info = 12;
  जइ (args.ldb < MAX(1,nrowa))  info =  9;
  जइ (args.lda < MAX(1,nrowa))  info =  7;
  जइ (args.k < 0)               info =  4;
  जइ (args.n < 0)               info =  3;
  जइ (trans < 0)                info =  2;
  जइ (uplo  < 0)                info =  1;

  जइ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अण्णहा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, गणना CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#अण्णहा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
#endif
#जइ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अण्णहा
	   शून्य *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#अण्णहा
	   शून्य *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#जइ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  पूर्णाङ्क uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अण्णहा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अण्णहा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अण्णहा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
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
#अण्णहा
  args.alpha = (शून्य *)alpha;
#endif

#जइ !defined(COMPLEX) || defined(HEMM)
  args.beta  = (शून्य *)&beta;
#अण्णहा
  args.beta  = (शून्य *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  जइ (order == CblasColMajor) {
    जइ (Uplo == CblasUpper) uplo  = 0;
    जइ (Uplo == CblasLower) uplo  = 1;

    जइ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    जइ (Trans == CblasTrans)       trans = 1;
    जइ (Trans == CblasConjNoTrans) trans = 0;
    जइ (Trans == CblasConjTrans)   trans = 1;
#नोचेद्यदि !defined(HEMM)
    जइ (Trans == CblasTrans)       trans = 1;
#अण्णहा
    जइ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    जइ (trans & 1) nrowa = args.k;

    जइ (args.ldc < MAX(1,args.n)) info = 12;
    जइ (args.ldb < MAX(1,nrowa))  info =  9;
    जइ (args.lda < MAX(1,nrowa))  info =  7;
    जइ (args.k < 0)               info =  4;
    जइ (args.n < 0)               info =  3;
    जइ (trans < 0)                info =  2;
    जइ (uplo  < 0)                info =  1;
  }

  जइ (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (शून्य *)CAlpha;
#endif

    जइ (Uplo == CblasUpper) uplo  = 1;
    जइ (Uplo == CblasLower) uplo  = 0;

    जइ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    जइ (Trans == CblasTrans)       trans = 0;
    जइ (Trans == CblasConjNoTrans) trans = 1;
    जइ (Trans == CblasConjTrans)   trans = 0;
#नोचेद्यदि !defined(HEMM)
    जइ (Trans == CblasTrans)       trans = 0;
#अण्णहा
    जइ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    जइ (trans & 1) nrowa = args.k;

    जइ (args.ldc < MAX(1,args.n)) info = 12;
    जइ (args.ldb < MAX(1,nrowa))  info =  9;
    जइ (args.lda < MAX(1,nrowa))  info =  7;
    जइ (args.k < 0)               info =  4;
    जइ (args.n < 0)               info =  3;
    जइ (trans < 0)                info =  2;
    जइ (uplo  < 0)                info =  1;
  }

  जइ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  जइ (args.n == 0) प्रत्यागम;

#ifdef DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#endif

#जइ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#जइ defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जइ defined(DYNAMIC_ARCH)
जइ (strcmp(gotoblas_corename(), "armv9sme") == 0
#जइ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   जइ (order == CblasRowMajor && n == ldc) {
     जइ (Trans == CblasNoTrans && k == lda && k == ldb) {
      जइ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अण्णहा जइ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }
     }
     अण्णहा जइ (Trans == CblasTrans && n == lda && n ==ldb) {
      जइ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अण्णहा जइ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  जइ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  जइ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अण्णहा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = रिक्त;
  जइ (args.n*args.k <1000)
	  args.nthreads =1 ;
  अण्णहा
  args.nthreads = अङ्कः_cpu_avail(3);

  जइ (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP

  } अण्णहा {

    syrk_thread(mode, &args, रिक्त, रिक्त, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  प्रत्यागम;
}
