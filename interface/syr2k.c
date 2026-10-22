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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <ctype.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QSYR2K"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DSYR2K"
#નહીંતર
#define ERROR_NAME "SSYR2K"
#endif
#નહીંતર
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYR2K"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZSYR2K"
#નહીંતર
#define ERROR_NAME "CSYR2K"
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XHER2K"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHER2K"
#નહીંતર
#define ERROR_NAME "CHER2K"
#endif
#endif
#endif

સ્થિર પૂર્ણાંક (*syr2k[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYR2K_UN, SYR2K_UC, SYR2K_LN, SYR2K_LC,
#નહીંતર
  HER2K_UN, HER2K_UC, HER2K_LN, HER2K_LC,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, અક્ષર *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
	               FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  અક્ષર uplo_arg  = *UPLO;
  અક્ષર trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  પૂર્ણાંક uplo;
  પૂર્ણાંક trans;
  પૂર્ણાંક nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (શૂન્ય *)a;
  args.b = (શૂન્ય *)b;
  args.c = (શૂન્ય *)c;

  args.lda = *ldA;
  args.ldb = *ldB;
  args.ldc = *ldC;

  args.alpha = (શૂન્ય *)alpha;
  args.beta  = (શૂન્ય *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;

#ifndef COMPLEX
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'T') trans = 1;
  જો (trans_arg == 'C') trans = 1;
#નહીંતર
#ifdef HEMM
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'C') trans = 1;
#નહીંતર
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'T') trans = 1;
#endif

#endif


  nrowa = args.n;
  જો (trans & 1) nrowa = args.k;

  info = 0;

  જો (args.ldc < MAX(1,args.n)) info = 12;
  જો (args.ldb < MAX(1,nrowa))  info =  9;
  જો (args.lda < MAX(1,nrowa))  info =  7;
  જો (args.k < 0)               info =  4;
  જો (args.n < 0)               info =  3;
  જો (trans < 0)                info =  2;
  જો (uplo  < 0)                info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_UPLO Uplo, ગણના CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
#નહીંતર
	   શૂન્ય *valpha,
	   શૂન્ય *va, blasint lda,
	   શૂન્ય *vb, blasint ldb,
#endif
#જો !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#નહીંતર
	   શૂન્ય *vbeta,
#endif
#ifndef COMPLEX
	   FLOAT *c,
#નહીંતર
	   શૂન્ય *vc,
#endif
           blasint ldc) {

#ifdef COMPLEX
  FLOAT* alpha = (FLOAT*) valpha;
#જો !defined(HEMM)
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* b = (FLOAT*) vb;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  પૂર્ણાંક uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef HEMM
  FLOAT CAlpha[2];
#endif

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (શૂન્ય *)a;
  args.b = (શૂન્ય *)b;
  args.c = (શૂન્ય *)c;

  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

#ifndef COMPLEX
  args.alpha = (શૂન્ય *)&alpha;
#નહીંતર
  args.alpha = (શૂન્ય *)alpha;
#endif

#જો !defined(COMPLEX) || defined(HEMM)
  args.beta  = (શૂન્ય *)&beta;
#નહીંતર
  args.beta  = (શૂન્ય *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (Uplo == CblasUpper) uplo  = 0;
    જો (Uplo == CblasLower) uplo  = 1;

    જો (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    જો (Trans == CblasTrans)       trans = 1;
    જો (Trans == CblasConjNoTrans) trans = 0;
    જો (Trans == CblasConjTrans)   trans = 1;
#અથવાજો !defined(HEMM)
    જો (Trans == CblasTrans)       trans = 1;
#નહીંતર
    જો (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    જો (trans & 1) nrowa = args.k;

    જો (args.ldc < MAX(1,args.n)) info = 12;
    જો (args.ldb < MAX(1,nrowa))  info =  9;
    જો (args.lda < MAX(1,nrowa))  info =  7;
    જો (args.k < 0)               info =  4;
    જો (args.n < 0)               info =  3;
    જો (trans < 0)                info =  2;
    જો (uplo  < 0)                info =  1;
  }

  જો (order == CblasRowMajor) {

#ifdef HEMM
    CAlpha[0] =  alpha[0];
    CAlpha[1] = -alpha[1];

    args.alpha = (શૂન્ય *)CAlpha;
#endif

    જો (Uplo == CblasUpper) uplo  = 1;
    જો (Uplo == CblasLower) uplo  = 0;

    જો (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    જો (Trans == CblasTrans)       trans = 0;
    જો (Trans == CblasConjNoTrans) trans = 1;
    જો (Trans == CblasConjTrans)   trans = 0;
#અથવાજો !defined(HEMM)
    જો (Trans == CblasTrans)       trans = 0;
#નહીંતર
    જો (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    જો (trans & 1) nrowa = args.k;

    જો (args.ldc < MAX(1,args.n)) info = 12;
    જો (args.ldb < MAX(1,nrowa))  info =  9;
    જો (args.lda < MAX(1,nrowa))  info =  7;
    જો (args.k < 0)               info =  4;
    જો (args.n < 0)               info =  3;
    જો (trans < 0)                info =  2;
    જો (uplo  < 0)                info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  જો (args.n == 0) પાછા;

#ifdef DYNAMIC_ARCH
બાહ્ય અક્ષર* gotoblas_corename(શૂન્ય);
#endif

#જો !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#જો defined(ARCH_ARM64) && (defined(USE_SSYR2K_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#જો defined(DYNAMIC_ARCH)
જો (strcmp(gotoblas_corename(), "armv9sme") == 0
#જો defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
   જો (order == CblasRowMajor && n == ldc) {
     જો (Trans == CblasNoTrans && k == lda && k == ldb) {
      જો (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        પાછા;
      }નહીંતર જો (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LN(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        પાછા;
      }
     }
     નહીંતર જો (Trans == CblasTrans && n == lda && n ==ldb) {
      જો (Uplo == CblasUpper) {
        SSYR2K_DIRECT_ALPHA_BETA_UT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        પાછા;
      }નહીંતર જો (Uplo == CblasLower) {
        SSYR2K_DIRECT_ALPHA_BETA_LT(n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        પાછા;
      }
     }
   }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  જો (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } નહીંતર {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }

  mode |= (uplo  << BLAS_UPLO_SHIFT);

  args.common = ખાલી;
  જો (args.n*args.k <1000)
	  args.nthreads =1 ;
  નહીંતર
  args.nthreads = અંક_cpu_avail(3);

  જો (args.nthreads == 1) {
#endif

    (syr2k[(uplo << 1) | trans ])(&args, ખાલી, ખાલી, sa, sb, 0);

#ifdef SMP

  } નહીંતર {

    syrk_thread(mode, &args, ખાલી, ખાલી, syr2k[(uplo << 1) | trans ], sa, sb, args.nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, 2 * args.n * args.k + args.n * args.n, 2 * args.n * args.n * args.k);

  IDEBUG_END;

  પાછા;
}
