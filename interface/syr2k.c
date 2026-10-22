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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#وإلا
#define ERROR_NAME "SSYR2K"
#endif
#وإلا
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#وإلا
#define ERROR_NAME "CSYR2K"
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#وإلا
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

ساكن صحيح (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#وإلا
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وإلا
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (فراغ *)a;
  args.b = (فراغ *)b;
  args.c = (فراغ *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  𐤀𐤌 (uplo_arg  == 'U') uplo  = 0;
  𐤀𐤌 (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  𐤀𐤌 (trans_arg == 'N') trans = 0;
  𐤀𐤌 (trans_arg == 'T') trans = 1;
  𐤀𐤌 (trans_arg == 'C') trans = 1;
#وإلا
#ifdef HEMM
  𐤀𐤌 (trans_arg == 'N') trans = 0;
  𐤀𐤌 (trans_arg == 'C') trans = 1;
#وإلا
  𐤀𐤌 (trans_arg == 'N') trans = 0;
  𐤀𐤌 (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  𐤀𐤌 (trans & 1) nrowa = args.k;

  info = 0;

  𐤀𐤌 (args.ldc < MAX(1,args.n)) info = 12;
  𐤀𐤌 (args.ldb < MAX(1,nrowa))  info =  9;
  𐤀𐤌 (args.lda < MAX(1,nrowa))  info =  7;
  𐤀𐤌 (args.k < 0)               info =  4;
  𐤀𐤌 (args.n < 0)               info =  3;
  𐤀𐤌 (trans < 0)                info =  2;
  𐤀𐤌 (uplo  < 0)                info =  1;

  𐤀𐤌 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
#endif
#𐤀𐤌 !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#وإلا
	   فراغ *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#𐤀𐤌 !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وإلا
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (فراغ *)a;
  args.b = (فراغ *)b;
  args.c = (فراغ *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (فراغ *)&alpha;
#وإلا
  args.alpha = (فراغ *)alpha;
#endif

#𐤀𐤌 !defined(COMPLEX) || defined(HEMM)
  args.beta  = (فراغ *)&beta;
#وإلا
  args.beta  = (فراغ *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  𐤀𐤌 (order == CblasColMajor) {
    𐤀𐤌 (Uplo == CblasUpper) uplo  = 0;
    𐤀𐤌 (Uplo == CblasLower) uplo  = 1;

    𐤀𐤌 (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    𐤀𐤌 (Trans == CblasTrans)       trans = 1;
    𐤀𐤌 (Trans == CblasConjNoTrans) trans = 0;
    𐤀𐤌 (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    𐤀𐤌 (Trans == CblasTrans)       trans = 1;
#وإلا
    𐤀𐤌 (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    𐤀𐤌 (trans & 1) nrowa = args.k;

    𐤀𐤌 (args.ldc < MAX(1,args.n)) info = 12;
    𐤀𐤌 (args.ldb < MAX(1,nrowa))  info =  9;
    𐤀𐤌 (args.lda < MAX(1,nrowa))  info =  7;
    𐤀𐤌 (args.k < 0)               info =  4;
    𐤀𐤌 (args.n < 0)               info =  3;
    𐤀𐤌 (trans < 0)                info =  2;
    𐤀𐤌 (uplo  < 0)                info =  1;
  }

  𐤀𐤌 (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (فراغ *)CAlpha;
#endif

    𐤀𐤌 (Uplo == CblasUpper) uplo  = 1;
    𐤀𐤌 (Uplo == CblasLower) uplo  = 0;

    𐤀𐤌 (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    𐤀𐤌 (Trans == CblasTrans)       trans = 0;
    𐤀𐤌 (Trans == CblasConjNoTrans) trans = 1;
    𐤀𐤌 (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    𐤀𐤌 (Trans == CblasTrans)       trans = 0;
#وإلا
    𐤀𐤌 (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    𐤀𐤌 (trans & 1) nrowa = args.k;

    𐤀𐤌 (args.ldc < MAX(1,args.n)) info = 12;
    𐤀𐤌 (args.ldb < MAX(1,nrowa))  info =  9;
    𐤀𐤌 (args.lda < MAX(1,nrowa))  info =  7;
    𐤀𐤌 (args.k < 0)               info =  4;
    𐤀𐤌 (args.n < 0)               info =  3;
    𐤀𐤌 (trans < 0)                info =  2;
    𐤀𐤌 (uplo  < 0)                info =  1;
  }

  𐤀𐤌 (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  𐤀𐤌 (args.n == 0) أرجع;

#ifdef DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#endif

#𐤀𐤌 !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#𐤀𐤌 defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#𐤀𐤌 defined(DYNAMIC_ARCH)
𐤀𐤌 (strcmp(gotoblas_corename(), "armv9sme") == 0
#𐤀𐤌 defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   𐤀𐤌 (order == CblasRowMajor && n == ldc) {
     𐤀𐤌 (Trans == CblasNoTrans && k == lda && k == ldb) {
      𐤀𐤌 (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا 𐤀𐤌 (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
     وإلا 𐤀𐤌 (Trans == CblasTrans && n == lda && n ==ldb) {
      𐤀𐤌 (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }وإلا 𐤀𐤌 (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        أرجع;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  𐤀𐤌 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  𐤀𐤌 (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = عدم;
  𐤀𐤌 (args.n*args.k <1000)
	  args.nthreads =1 ;
  وإلا
  args.nthreads = num_cpu_avail(3);

  𐤀𐤌 (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP

  } وإلا {

    syrk_thread(mode, &args, عدم, عدم, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  أرجع;
}
