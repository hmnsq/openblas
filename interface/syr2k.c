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
#अन्यथा
#define ERROR_NAME "SSYR2K"
#endif
#अन्यथा
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#अन्यथा
#define ERROR_NAME "CSYR2K"
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#अन्यथा
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

स्थिर पूर्णाङ्क (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#अन्यथा
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
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
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

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
#ifdef HEMM
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  यदि (trans & 1) nrowa = args.k;

  info = 0;

  यदि (args.ldc < MAX(1,args.n)) info = 12;
  यदि (args.ldb < MAX(1,nrowa))  info =  9;
  यदि (args.lda < MAX(1,nrowa))  info =  7;
  यदि (args.k < 0)               info =  4;
  यदि (args.n < 0)               info =  3;
  यदि (trans < 0)                info =  2;
  यदि (uplo  < 0)                info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
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
#यदि !defined(COMPLEX) || defined(HEMM)
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
#यदि !defined(HEMM)
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
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
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
#अन्यथा
  args.alpha = (शून्य *)alpha;
#endif

#यदि !defined(COMPLEX) || defined(HEMM)
  args.beta  = (शून्य *)&beta;
#अन्यथा
  args.beta  = (शून्य *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  यदि (order == CblasColMajor) {
    यदि (Uplo == CblasUpper) uplo  = 0;
    यदि (Uplo == CblasLower) uplo  = 1;

    यदि (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    यदि (Trans == CblasTrans)       trans = 1;
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 1;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    यदि (trans & 1) nrowa = args.k;

    यदि (args.ldc < MAX(1,args.n)) info = 12;
    यदि (args.ldb < MAX(1,nrowa))  info =  9;
    यदि (args.lda < MAX(1,nrowa))  info =  7;
    यदि (args.k < 0)               info =  4;
    यदि (args.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (शून्य *)CAlpha;
#endif

    यदि (Uplo == CblasUpper) uplo  = 1;
    यदि (Uplo == CblasLower) uplo  = 0;

    यदि (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    यदि (Trans == CblasTrans)       trans = 0;
    यदि (Trans == CblasConjNoTrans) trans = 1;
    यदि (Trans == CblasConjTrans)   trans = 0;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 0;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    यदि (trans & 1) nrowa = args.k;

    यदि (args.ldc < MAX(1,args.n)) info = 12;
    यदि (args.ldb < MAX(1,nrowa))  info =  9;
    यदि (args.lda < MAX(1,nrowa))  info =  7;
    यदि (args.k < 0)               info =  4;
    यदि (args.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (args.n == 0) प्रत्यागम;

#ifdef DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#endif

#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (strcmp(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   यदि (order == CblasRowMajor && n == ldc) {
     यदि (Trans == CblasNoTrans && k == lda && k == ldb) {
      यदि (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अन्यथा यदि (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }
     }
     अन्यथा यदि (Trans == CblasTrans && n == lda && n ==ldb) {
      यदि (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        प्रत्यागम;
      }अन्यथा यदि (Uplo == CblasLower) {
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
  यदि (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  यदि (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = रिक्त;
  यदि (args.n*args.k <1000)
	  args.nthreads =1 ;
  अन्यथा
  args.nthreads = अंक_cpu_avail(3);

  यदि (args.nthreads == 1) {
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

  प्रत्यागम;
}
