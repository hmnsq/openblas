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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#elif defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#ပေႃးဢမ်ႇ
#define ERROR_NAME "STRSM "
#endif
#ပေႃးဢမ်ႇ
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#ပေႃးဢမ်ႇ
#define ERROR_NAME "CTRSM "
#endif
#endif
#ပေႃးဢမ်ႇ
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#elif defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#ပေႃးဢမ်ႇ
#define ERROR_NAME "STRMM "
#endif
#ပေႃးဢမ်ႇ
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#ပေႃးဢမ်ႇ
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#ပေႃးဢမ်ႇ
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
extern var_na* gotoblas_corename(void);                                
#endif                                                               

static puur_naa_nka (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#ပေႃးဢမ်ႇ
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *SIDE, var_na *UPLO, var_na *TRANS, var_na *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  var_na side_arg  = *SIDE;
  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;
  var_na diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  puur_naa_nka mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  puur_naa_nka mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ပေႃးဢမ်ႇ
  puur_naa_nka mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ပေႃးဢမ်ႇ
#ifdef XDOUBLE
  puur_naa_nka mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  puur_naa_nka mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ပေႃးဢမ်ႇ
  puur_naa_nka mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  puur_naa_nka side;
  puur_naa_nka uplo;
  puur_naa_nka unit;
  puur_naa_nka trans;
  puur_naa_nka nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (void *)a;
  args.b = (void *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (void *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  ပေႃး (side_arg  == 'L') side  = 0;
  ပေႃး (side_arg  == 'R') side  = 1;

  ပေႃး (trans_arg == 'N') trans = 0;
  ပေႃး (trans_arg == 'T') trans = 1;
  ပေႃး (trans_arg == 'R') trans = 2;
  ပေႃး (trans_arg == 'C') trans = 3;

  ပေႃး (diag_arg  == 'U') unit  = 0;
  ပေႃး (diag_arg  == 'N') unit  = 1;

  ပေႃး (uplo_arg  == 'U') uplo  = 0;
  ပေႃး (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  ပေႃး (side & 1) nrowa = args.n;

  info = 0;

  ပေႃး (args.ldb < MAX(1,args.m)) info = 11;
  ပေႃး (args.lda < MAX(1,nrowa))  info =  9;
  ပေႃး (args.n < 0)               info =  6;
  ပေႃး (args.m < 0)               info =  5;
  ပေႃး (unit < 0)                 info =  4;
  ပေႃး (trans < 0)                info =  3;
  ပေႃး (uplo  < 0)                info =  2;
  ပေႃး (side  < 0)                info =  1;

  ပေႃး (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME)-1);
    ပွၵ်ႈ;
  }

#ပေႃးဢမ်ႇ

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_SIDE Side,  ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE Trans, ga_nanaa CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#ပေႃးဢမ်ႇ
	   void *valpha,
	   void *va, blasint lda,
	   void *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#endif

  blas_arg_t args;
  puur_naa_nka side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  puur_naa_nka mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  puur_naa_nka mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ပေႃးဢမ်ႇ
  puur_naa_nka mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ပေႃးဢမ်ႇ
#ifdef XDOUBLE
  puur_naa_nka mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  puur_naa_nka mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ပေႃးဢမ်ႇ
  puur_naa_nka mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (void *)a;
  args.b = (void *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (void *)&alpha;
#ပေႃးဢမ်ႇ
  args.beta = (void *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  ပေႃး (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    ပေႃး (Side == CblasLeft)         side  = 0;
    ပေႃး (Side == CblasRight)        side  = 1;

    ပေႃး (Uplo == CblasUpper)        uplo  = 0;
    ပေႃး (Uplo == CblasLower)        uplo  = 1;

    ပေႃး (Trans == CblasNoTrans)     trans = 0;
    ပေႃး (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    ပေႃး (Trans == CblasConjNoTrans) trans = 0;
    ပေႃး (Trans == CblasConjTrans)   trans = 1;
#ပေႃးဢမ်ႇ
    ပေႃး (Trans == CblasConjNoTrans) trans = 2;
    ပေႃး (Trans == CblasConjTrans)   trans = 3;
#endif

    ပေႃး (Diag == CblasUnit)          unit  = 0;
    ပေႃး (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    ပေႃး (side & 1) nrowa = args.n;

    ပေႃး (args.ldb < MAX(1,args.m)) info = 11;
    ပေႃး (args.lda < MAX(1,nrowa))  info =  9;
    ပေႃး (args.n < 0)               info =  6;
    ပေႃး (args.m < 0)               info =  5;
    ပေႃး (unit < 0)                 info =  4;
    ပေႃး (trans < 0)                info =  3;
    ပေႃး (uplo  < 0)                info =  2;
    ပေႃး (side  < 0)                info =  1;
  }

  ပေႃး (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    ပေႃး (Side == CblasLeft)         side  = 1;
    ပေႃး (Side == CblasRight)        side  = 0;

    ပေႃး (Uplo == CblasUpper)        uplo  = 1;
    ပေႃး (Uplo == CblasLower)        uplo  = 0;

    ပေႃး (Trans == CblasNoTrans)     trans = 0;
    ပေႃး (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    ပေႃး (Trans == CblasConjNoTrans) trans = 0;
    ပေႃး (Trans == CblasConjTrans)   trans = 1;
#ပေႃးဢမ်ႇ
    ပေႃး (Trans == CblasConjNoTrans) trans = 2;
    ပေႃး (Trans == CblasConjTrans)   trans = 3;
#endif

    ပေႃး (Diag == CblasUnit)         unit  = 0;
    ပေႃး (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    ပေႃး (side & 1) nrowa = args.n;

    ပေႃး (args.ldb < MAX(1,args.m)) info = 11;
    ပေႃး (args.lda < MAX(1,nrowa))  info =  9;
    ပေႃး (args.n < 0)               info =  6;
    ပေႃး (args.m < 0)               info =  5;
    ပေႃး (unit < 0)                 info =  4;
    ပေႃး (trans < 0)                info =  3;
    ပေႃး (uplo  < 0)                info =  2;
    ပေႃး (side  < 0)                info =  1;
  }

  ပေႃး (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

  ပေႃး (args.m == 0 || args.n == 0) ပွၵ်ႈ;

#ပေႃး defined(TRMM)
#ပေႃး !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ပေႃး defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ပေႃး defined(DYNAMIC_ARCH)
ပေႃး (strcmp(gotoblas_corename(), "armv9sme") == 0
#ပေႃး defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  ပေႃး (order == CblasRowMajor && Diag == CblasNonUnit && Side == CblasLeft && m == lda && n == ldb) {
    ပေႃး (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      ပွၵ်ႈ;
    } ပေႃးဢမ်ႇ ပေႃး (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      ပွၵ်ႈ;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  ပေႃး ((args.m == 0) || (args.n == 0)) ပွၵ်ႈ;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ပေႃး (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပွၵ်ႈ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  mode |= (trans << BLAS_TRANSA_SHIFT);
  mode |= (side  << BLAS_RSIDE_SHIFT);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  ပေႃး ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  ပေႃးဢမ်ႇ
	args.nthreads = num_cpu_avail(3);
		

  ပေႃး (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } ပေႃးဢမ်ႇ {
    ပေႃး (!side) {
      gemm_thread_n(mode, &args, NULL, NULL, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } ပေႃးဢမ်ႇ {
      gemm_thread_m(mode, &args, NULL, NULL, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  ပွၵ်ႈ;
}

