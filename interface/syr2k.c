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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <ctype.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#അല്ലെങ്കിൽ
#define ERROR_NAME "SSYR2K"
#endif
#അല്ലെങ്കിൽ
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#അല്ലെങ്കിൽ
#define ERROR_NAME "CSYR2K"
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#അല്ലെങ്കിൽ
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#അല്ലെങ്കിൽ
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(അക്ഷരം *UPLO, അക്ഷരം *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  അക്ഷരം uplo_arg  = *UPLO;
  അക്ഷരം trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  പൂർണ്ണസംഖ്യ uplo;
  പൂർണ്ണസംഖ്യ trans;
  പൂർണ്ണസംഖ്യ nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (ശൂന്യം *)a;
  args.b = (ശൂന്യം *)b;
  args.c = (ശൂന്യം *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (ശൂന്യം *)alpha;
  args.beta  = (ശൂന്യം *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  എങ്കിൽ (uplo_arg  == 'U') uplo  = 0;
  എങ്കിൽ (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'T') trans = 1;
  എങ്കിൽ (trans_arg == 'C') trans = 1;
#അല്ലെങ്കിൽ
#ifdef HEMM
  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'C') trans = 1;
#അല്ലെങ്കിൽ
  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  എങ്കിൽ (trans & 1) nrowa = args.k;

  info = 0;

  എങ്കിൽ (args.ldc < MAX(1,args.n)) info = 12;
  എങ്കിൽ (args.ldb < MAX(1,nrowa))  info =  9;
  എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  7;
  എങ്കിൽ (args.k < 0)               info =  4;
  എങ്കിൽ (args.n < 0)               info =  3;
  എങ്കിൽ (trans < 0)                info =  2;
  എങ്കിൽ (uplo  < 0)                info =  1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_UPLO Uplo, എണ്ണൽ CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#അല്ലെങ്കിൽ
	   ശൂന്യം *valpha,
	   ശൂന്യം *va, blasint lda,
	   ശൂന്യം *vb, blasint ldb,
#endif
#എങ്കിൽ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#അല്ലെങ്കിൽ
	   ശൂന്യം *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#അല്ലെങ്കിൽ
	   ശൂന്യം *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#എങ്കിൽ !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  പൂർണ്ണസംഖ്യ uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  പൂർണ്ണസംഖ്യ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  പൂർണ്ണസംഖ്യ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  പൂർണ്ണസംഖ്യ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (ശൂന്യം *)a;
  args.b = (ശൂന്യം *)b;
  args.c = (ശൂന്യം *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (ശൂന്യം *)&alpha;
#അല്ലെങ്കിൽ
  args.alpha = (ശൂന്യം *)alpha;
#endif

#എങ്കിൽ !defined(COMPLEX) || defined(HEMM)
  args.beta  = (ശൂന്യം *)&beta;
#അല്ലെങ്കിൽ
  args.beta  = (ശൂന്യം *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  എങ്കിൽ (order == CblasColMajor) {
    എങ്കിൽ (Uplo == CblasUpper) uplo  = 0;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 1;

    എങ്കിൽ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    എങ്കിൽ (Trans == CblasTrans)       trans = 1;
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 0;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 1;
#അല്ലെങ്കിലെങ്കിൽ !defined(HEMM)
    എങ്കിൽ (Trans == CblasTrans)       trans = 1;
#അല്ലെങ്കിൽ
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    എങ്കിൽ (trans & 1) nrowa = args.k;

    എങ്കിൽ (args.ldc < MAX(1,args.n)) info = 12;
    എങ്കിൽ (args.ldb < MAX(1,nrowa))  info =  9;
    എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  7;
    എങ്കിൽ (args.k < 0)               info =  4;
    എങ്കിൽ (args.n < 0)               info =  3;
    എങ്കിൽ (trans < 0)                info =  2;
    എങ്കിൽ (uplo  < 0)                info =  1;
  }

  എങ്കിൽ (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (ശൂന്യം *)CAlpha;
#endif

    എങ്കിൽ (Uplo == CblasUpper) uplo  = 1;
    എങ്കിൽ (Uplo == CblasLower) uplo  = 0;

    എങ്കിൽ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    എങ്കിൽ (Trans == CblasTrans)       trans = 0;
    എങ്കിൽ (Trans == CblasConjNoTrans) trans = 1;
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 0;
#അല്ലെങ്കിലെങ്കിൽ !defined(HEMM)
    എങ്കിൽ (Trans == CblasTrans)       trans = 0;
#അല്ലെങ്കിൽ
    എങ്കിൽ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    എങ്കിൽ (trans & 1) nrowa = args.k;

    എങ്കിൽ (args.ldc < MAX(1,args.n)) info = 12;
    എങ്കിൽ (args.ldb < MAX(1,nrowa))  info =  9;
    എങ്കിൽ (args.lda < MAX(1,nrowa))  info =  7;
    എങ്കിൽ (args.k < 0)               info =  4;
    എങ്കിൽ (args.n < 0)               info =  3;
    എങ്കിൽ (trans < 0)                info =  2;
    എങ്കിൽ (uplo  < 0)                info =  1;
  }

  എങ്കിൽ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

  എങ്കിൽ (args.n == 0) തിരികെ;

#ifdef DYNAMIC_ARCH
ബാഹ്യ അക്ഷരം* gotoblas_corename(ശൂന്യം);
#endif

#എങ്കിൽ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#എങ്കിൽ defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#എങ്കിൽ defined(DYNAMIC_ARCH)
എങ്കിൽ (strcmp(gotoblas_corename(), "armv9sme") == 0
#എങ്കിൽ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   എങ്കിൽ (order == CblasRowMajor && n == ldc) {
     എങ്കിൽ (Trans == CblasNoTrans && k == lda && k == ldb) {
      എങ്കിൽ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        തിരികെ;
      }അല്ലെങ്കിൽ എങ്കിൽ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        തിരികെ;
      }
     }
     അല്ലെങ്കിൽ എങ്കിൽ (Trans == CblasTrans && n == lda && n ==ldb) {
      എങ്കിൽ (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        തിരികെ;
      }അല്ലെങ്കിൽ എങ്കിൽ (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        തിരികെ;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  എങ്കിൽ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
    തിരികെ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  എങ്കിൽ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } അല്ലെങ്കിൽ {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = ശൂന്യ;
  എങ്കിൽ (args.n*args.k <1000)
	  args.nthreads =1 ;
  അല്ലെങ്കിൽ
  args.nthreads = അക്കം_cpu_avail(3);

  എങ്കിൽ (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP

  } അല്ലെങ്കിൽ {

    syrk_thread(mode, &args, ശൂന്യ, ശൂന്യ, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  തിരികെ;
}
