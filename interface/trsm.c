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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_نه_تعریف_شوی TRMM
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DTRSM "
#وإلا
#تعریف_کړه ERROR_NAME "STRSM "
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XTRSM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZTRSM "
#وإلا
#تعریف_کړه ERROR_NAME "CTRSM "
#د_که_پای
#د_که_پای
#وإلا
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DTRMM "
#وإلا
#تعریف_کړه ERROR_NAME "STRMM "
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XTRMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZTRMM "
#وإلا
#تعریف_کړه ERROR_NAME "CTRMM "
#د_که_پای
#د_که_پای
#د_که_پای

#که_نه_تعریف_شوی COMPLEX
#تعریف_کړه SMP_ضریب 256
#وإلا
#تعریف_کړه SMP_ضریب 128
#د_که_پای

#که_تعریف_شوی DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#د_که_پای                                                               

ساكن صحيح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#که_نه_تعریف_شوی TRMM
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
#د_که_پای
};

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blas_arg_t دلیلونه;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای

  blasint info;
  صحيح side;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دلیلونه.m = *M;
  دلیلونه.n = *N;

  دلیلونه.a = (فراغ *)a;
  دلیلونه.b = (فراغ *)b;

  دلیلونه.lda = *ldA;
  دلیلونه.ldb = *ldB;

  دلیلونه.beta = (فراغ *)alpha;

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

  nrowa = دلیلونه.m;
  إذا (side & 1) nrowa = دلیلونه.n;

  info = 0;

  إذا (دلیلونه.ldb < MAX(1,دلیلونه.m)) info = 11;
  إذا (دلیلونه.lda < MAX(1,nrowa))  info =  9;
  إذا (دلیلونه.n < 0)               info =  6;
  إذا (دلیلونه.m < 0)               info =  5;
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
#که_نه_تعریف_شوی COMPLEX
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
#د_که_پای

  blas_arg_t دلیلونه;
  صحيح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  صحيح ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای

  PRINT_DEBUG_CNAME;

  دلیلونه.a = (فراغ *)a;
  دلیلونه.b = (فراغ *)b;

  دلیلونه.lda = lda;
  دلیلونه.ldb = ldb;

#که_نه_تعریف_شوی COMPLEX
  دلیلونه.beta = (فراغ *)&alpha;
#وإلا
  دلیلونه.beta = (فراغ *)alpha;
#د_که_پای

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    دلیلونه.m = m;
    دلیلونه.n = n;

    إذا (Side == Cblas_کیڼ)         side  = 0;
    إذا (Side == Cblas_ښی)        side  = 1;

    إذا (Uplo == Cblas_پاسنی)        uplo  = 0;
    إذا (Uplo == Cblas_ښکتنی)        uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#د_که_پای

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = دلیلونه.m;
    إذا (side & 1) nrowa = دلیلونه.n;

    إذا (دلیلونه.ldb < MAX(1,دلیلونه.m)) info = 11;
    إذا (دلیلونه.lda < MAX(1,nrowa))  info =  9;
    إذا (دلیلونه.n < 0)               info =  6;
    إذا (دلیلونه.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    دلیلونه.m = n;
    دلیلونه.n = m;

    إذا (Side == Cblas_کیڼ)         side  = 1;
    إذا (Side == Cblas_ښی)        side  = 0;

    إذا (Uplo == Cblas_پاسنی)        uplo  = 1;
    إذا (Uplo == Cblas_ښکتنی)        uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 0;
    إذا (Trans == CblasTrans)       trans = 1;
#که_نه_تعریف_شوی COMPLEX
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلا
    إذا (Trans == CblasConjNoTrans) trans = 2;
    إذا (Trans == CblasConjTrans)   trans = 3;
#د_که_پای

    إذا (Diag == CblasUnit)         unit  = 0;
    إذا (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = دلیلونه.m;
    إذا (side & 1) nrowa = دلیلونه.n;

    إذا (دلیلونه.ldb < MAX(1,دلیلونه.m)) info = 11;
    إذا (دلیلونه.lda < MAX(1,nrowa))  info =  9;
    إذا (دلیلونه.n < 0)               info =  6;
    إذا (دلیلونه.m < 0)               info =  5;
    إذا (unit < 0)                 info =  4;
    إذا (trans < 0)                info =  3;
    إذا (uplo  < 0)                info =  2;
    إذا (side  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دلیلونه.m == 0 || دلیلونه.n == 0) أرجع;

#إذا defined(TRMM)
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#د_که_پای
)
#د_که_پای
  إذا (ترتیب == Cblas_لیکه_Major && Diag == CblasNonUnit && Side == Cblas_کیڼ && m == lda && n == ldb) {
    إذا (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_پاسنی ? STRMM_مستقیم_LNUN : STRMM_مستقیم_LNLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    } وإلا إذا (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_پاسنی ? STRMM_مستقیم_LTUN : STRMM_مستقیم_LTLN)(m, n, alpha, a, lda, b, ldb);
      أرجع;
    }
  }
#د_که_پای
#د_که_پای
#د_که_پای



//end of the ifndef CBLAS ... else ... section	
#د_که_پای

  إذا ((دلیلونه.m == 0) || (دلیلونه.n == 0)) أرجع;
  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  بفر = (FLOAT *)blas_حافظه_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_دننه_واټن_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_دننه_واټن_B);

#که_تعریف_شوی SMP
  ډېرپېښ |= (trans << BLAS_TRANSA_ښويدنه);
  ډېرپېښ |= (side  << BLAS_RSIDE_ښويدنه);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  إذا ( دلیلونه.m * دلیلونه.n < SMP_ضریب * GEMM_MULTITHREAD_پوله)
	دلیلونه.nthreads = 1;
  وإلا
	دلیلونه.nthreads = عدد_cpu_avail(3);
		

  إذا (دلیلونه.nthreads == 1) {
#د_که_پای

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&دلیلونه, عدم, عدم, sa, sb, 0);

#که_تعریف_شوی SMP
  } وإلا {
    إذا (!side) {
      gemm_thread_n(ډېرپېښ, &دلیلونه, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلیلونه.nthreads);
    } وإلا {
      gemm_thread_m(ډېرپېښ, &دلیلونه, عدم, عدم, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلیلونه.nthreads);
    }
  }
#د_که_پای

  blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE * COMPSIZE,
		       (!side) ? دلیلونه.m * (دلیلونه.m + دلیلونه.n) : دلیلونه.n * (دلیلونه.m + دلیلونه.n),
		       (!side) ? دلیلونه.m * دلیلونه.m * دلیلونه.n : دلیلونه.m * دلیلونه.n * دلیلونه.n);

  IDEBUG_پای_ټکی;

  أرجع;
}

