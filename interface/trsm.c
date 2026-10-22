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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل TRMM
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QTRSM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DTRSM "
#وإلا
#وصف_ڪريو ERROR_NAME "STRSM "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XTRSM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZTRSM "
#وإلا
#وصف_ڪريو ERROR_NAME "CTRSM "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QTRMM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DTRMM "
#وإلا
#وصف_ڪريو ERROR_NAME "STRMM "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XTRMM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZTRMM "
#وإلا
#وصف_ڪريو ERROR_NAME "CTRMM "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#وصف_ڪريو SMP_ضرب_ڪندڙ 256
#وإلا
#وصف_ڪريو SMP_ضرب_ڪندڙ 128
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#جيڪڏهن_پورو                                                               

ساكن صحيح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جيڪڏهن_وصف_نه_ٿيل TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#وإلا
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#جيڪڏهن_پورو
};

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blas_arg_t دليلون;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  blasint info;
  صحيح side;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دليلون.m = *M;
  دليلون.n = *N;

  دليلون.a = (فراغ *)a;
  دليلون.b = (فراغ *)b;

  دليلون.lda = *ldA;
  دليلون.ldb = *ldB;

  دليلون.beta = (فراغ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (diag_arg  == 'U') unit  = 0;
  إذا (diag_arg  == 'N') unit  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  nrowa = دليلون.m;
  إذا (side & 1) nrowa = دليلون.n;

  info = 0;

  إذا (دليلون.ldb < MAX(1,دليلون.m)) info = 11;
  إذا (دليلون.lda < MAX(1,nrowa))  info =  9;
  إذا (دليلون.n < 0)               info =  6;
  إذا (دليلون.m < 0)               info =  5;
  إذا (unit < 0)                 info =  4;
  إذا (trans < 0)                info =  3;
  إذا (uplo  < 0)                info =  2;
  إذا (side  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME)-1);
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب,
	   تعداد CBLAS_SIDE Side,  تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE Trans, تعداد CBLAS_DIAG Diag,
	   blasint m, blasint n,
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#جيڪڏهن_پورو

  blas_arg_t دليلون;
  صحيح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  دليلون.a = (فراغ *)a;
  دليلون.b = (فراغ *)b;

  دليلون.lda = lda;
  دليلون.ldb = ldb;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  دليلون.beta = (فراغ *)&alpha;
#وإلا
  دليلون.beta = (فراغ *)alpha;
#جيڪڏهن_پورو

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    دليلون.m = m;
    دليلون.n = n;

    إذا (Side == Cblas_کاٻو)         side  = 0;
    إذا (Side == Cblas_ساڄو)        side  = 1;

    إذا (Uplo == Cblas_مٿانهون)        uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون)        uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#جيڪڏهن_پورو

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = دليلون.m;
    إذا (side & 1) nrowa = دليلون.n;

    إذا (دليلون.ldb < MAX(1,دليلون.m)) info = 11;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  9;
    إذا (دليلون.n < 0)               info =  6;
    إذا (دليلون.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
    دليلون.m = n;
    دليلون.n = m;

    إذا (Side == Cblas_کاٻو)         side  = 1;
    إذا (Side == Cblas_ساڄو)        side  = 0;

    إذا (Uplo == Cblas_مٿانهون)        uplo  = 1;
    إذا (Uplo == Cblas_هيٺاهون)        uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#جيڪڏهن_پورو

    إذا (Diag == CblasUnit)         unit  = 0;
    إذا (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = دليلون.m;
    إذا (side & 1) nrowa = دليلون.n;

    إذا (دليلون.ldb < MAX(1,دليلون.m)) info = 11;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  9;
    إذا (دليلون.n < 0)               info =  6;
    إذا (دليلون.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دليلون.m == 0 || دليلون.n == 0) أرجع;

#إذا defined(TRMM)
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_سڌو)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جيڪڏهن_پورو
)
#جيڪڏهن_پورو
  إذا (ترتيب == Cblas_صف_Major && Diag == CblasNonUnit && Side == Cblas_کاٻو && m == lda && n == ldb) {
    إذا (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_مٿانهون ? STRMM_سڌو_LNUN : STRMM_سڌو_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا إذا (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_مٿانهون ? STRMM_سڌو_LTUN : STRMM_سڌو_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو



//end of the ifndef CBLAS ... else ... section	
#جيڪڏهن_پورو

  إذا ((دليلون.m == 0) || (دليلون.n == 0)) أرجع;
  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  بفر = (FLOAT *)blas_يادگيري_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);

#جيڪڏهن_وصف_ٿيل SMP
  گهڻائي |= (trans << BLAS_TRANSA_سرڪايو);
  گهڻائي |= (side  << BLAS_RSIDE_سرڪايو);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  إذا ( دليلون.m * دليلون.n < SMP_ضرب_ڪندڙ * GEMM_MULTITHREAD_حد)
	دليلون.nthreads = 1;
  وإلا
	دليلون.nthreads = انگ_cpu_avail(3);
		

  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {
    إذا (!side) {
      gemm_thread_n(گهڻائي, &دليلون, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دليلون.nthreads);
    } وإلا {
      gemm_thread_m(گهڻائي, &دليلون, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دليلون.nthreads);
    }
  }
#جيڪڏهن_پورو

  blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE,
		       (!side) ? دليلون.m * (دليلون.m + دليلون.n) : دليلون.n * (دليلون.m + دليلون.n),
		       (!side) ? دليلون.m * دليلون.m * دليلون.n : دليلون.m * دليلون.n * دليلون.n);

  IDEBUG_پڇاڙي;

  أرجع;
}

