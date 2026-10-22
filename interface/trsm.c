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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_بنا_تعریف TRMM
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTRSM "
#وإلا
#تعریف_کرو ERROR_NAME "STRSM "
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZTRSM "
#وإلا
#تعریف_کرو ERROR_NAME "CTRSM "
#جے_انت
#جے_انت
#وإلا
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTRMM "
#وإلا
#تعریف_کرو ERROR_NAME "STRMM "
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZTRMM "
#وإلا
#تعریف_کرو ERROR_NAME "CTRMM "
#جے_انت
#جے_انت
#جے_انت

#جے_بنا_تعریف COMPLEX
#تعریف_کرو SMP_ضارب 256
#وإلا
#تعریف_کرو SMP_ضارب 128
#جے_انت

#جے_تعریف_شدہ DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#جے_انت                                                               

ساكن صحيح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جے_بنا_تعریف TRMM
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
#جے_انت
};

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blas_arg_t دلیلاں;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت

  blasint info;
  صحيح side;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دلیلاں.m = *M;
  دلیلاں.n = *N;

  دلیلاں.a = (فراغ *)a;
  دلیلاں.b = (فراغ *)b;

  دلیلاں.lda = *ldA;
  دلیلاں.ldb = *ldB;

  دلیلاں.beta = (فراغ *)alpha;

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

  nrowa = دلیلاں.m;
  إذا (side & 1) nrowa = دلیلاں.n;

  info = 0;

  إذا (دلیلاں.ldb < MAX(1,دلیلاں.m)) info = 11;
  إذا (دلیلاں.lda < MAX(1,nrowa))  info =  9;
  إذا (دلیلاں.n < 0)               info =  6;
  إذا (دلیلاں.m < 0)               info =  5;
  إذا (unit < 0)                 info =  4;
  إذا (trans < 0)                info =  3;
  إذا (uplo  < 0)                info =  2;
  إذا (side  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME)-1);
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب,
	   تعداد CBLAS_SIDE Side,  تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE Trans, تعداد CBLAS_DIAG Diag,
	   blasint m, blasint n,
#جے_بنا_تعریف COMPLEX
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
#جے_انت

  blas_arg_t دلیلاں;
  صحيح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت

  PRINT_DEBUG_CNAME;

  دلیلاں.a = (فراغ *)a;
  دلیلاں.b = (فراغ *)b;

  دلیلاں.lda = lda;
  دلیلاں.ldb = ldb;

#جے_بنا_تعریف COMPLEX
  دلیلاں.beta = (فراغ *)&alpha;
#وإلا
  دلیلاں.beta = (فراغ *)alpha;
#جے_انت

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {
    دلیلاں.m = m;
    دلیلاں.n = n;

    إذا (Side == Cblas_بایاں)         side  = 0;
    إذا (Side == Cblas_دایاں)        side  = 1;

    إذا (Uplo == Cblas_اتلا)        uplo  = 0;
    إذا (Uplo == Cblas_ہیٹھلا)        uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#جے_بنا_تعریف COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#جے_انت

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = دلیلاں.m;
    إذا (side & 1) nrowa = دلیلاں.n;

    إذا (دلیلاں.ldb < MAX(1,دلیلاں.m)) info = 11;
    إذا (دلیلاں.lda < MAX(1,nrowa))  info =  9;
    إذا (دلیلاں.n < 0)               info =  6;
    إذا (دلیلاں.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_صف_Major) {
    دلیلاں.m = n;
    دلیلاں.n = m;

    إذا (Side == Cblas_بایاں)         side  = 1;
    إذا (Side == Cblas_دایاں)        side  = 0;

    إذا (Uplo == Cblas_اتلا)        uplo  = 1;
    إذا (Uplo == Cblas_ہیٹھلا)        uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#جے_بنا_تعریف COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#جے_انت

    إذا (Diag == CblasUnit)         unit  = 0;
    إذا (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = دلیلاں.m;
    إذا (side & 1) nrowa = دلیلاں.n;

    إذا (دلیلاں.ldb < MAX(1,دلیلاں.m)) info = 11;
    إذا (دلیلاں.lda < MAX(1,nrowa))  info =  9;
    إذا (دلیلاں.n < 0)               info =  6;
    إذا (دلیلاں.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دلیلاں.m == 0 || دلیلاں.n == 0) أرجع;

#إذا defined(TRMM)
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جے_انت
)
#جے_انت
  إذا (ترتیب == Cblas_صف_Major && Diag == CblasNonUnit && Side == Cblas_بایاں && m == lda && n == ldb) {
    إذا (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_اتلا ? STRMM_سدھا_LNUN : STRMM_سدھا_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا إذا (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_اتلا ? STRMM_سدھا_LTUN : STRMM_سدھا_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#جے_انت
#جے_انت
#جے_انت



//end of the ifndef CBLAS ... else ... section	
#جے_انت

  إذا ((دلیلاں.m == 0) || (دلیلاں.n == 0)) أرجع;
  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وتھ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ناپ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وتھ_B);

#جے_تعریف_شدہ SMP
  منوال |= (trans << BLAS_TRANSA_کھسکاؤ);
  منوال |= (side  << BLAS_RSIDE_کھسکاؤ);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  إذا ( دلیلاں.m * دلیلاں.n < SMP_ضارب * GEMM_MULTITHREAD_حد)
	دلیلاں.nthreads = 1;
  وإلا
	دلیلاں.nthreads = عدد_cpu_avail(3);
		

  إذا (دلیلاں.nthreads == 1) {
#جے_انت

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&دلیلاں, عدم, عدم, sa, sb, 0);

#جے_تعریف_شدہ SMP
  } وإلا {
    إذا (!side) {
      gemm_thread_n(منوال, &دلیلاں, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلیلاں.nthreads);
    } وإلا {
      gemm_thread_m(منوال, &دلیلاں, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلیلاں.nthreads);
    }
  }
#جے_انت

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_سرا(COMPSIZE * COMPSIZE,
		       (!side) ? دلیلاں.m * (دلیلاں.m + دلیلاں.n) : دلیلاں.n * (دلیلاں.m + دلیلاں.n),
		       (!side) ? دلیلاں.m * دلیلاں.m * دلیلاں.n : دلیلاں.m * دلیلاں.n * دلیلاں.n);

  IDEBUG_سرا;

  أرجع;
}

