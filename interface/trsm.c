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

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_غیر_تعریف_شدہ TRMM
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTRSM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTRSM "
#ورنہ
#تعریف_کرو ERROR_NAME "STRSM "
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XTRSM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZTRSM "
#ورنہ
#تعریف_کرو ERROR_NAME "CTRSM "
#اگر_ختم
#اگر_ختم
#ورنہ
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTRMM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTRMM "
#ورنہ
#تعریف_کرو ERROR_NAME "STRMM "
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XTRMM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZTRMM "
#ورنہ
#تعریف_کرو ERROR_NAME "CTRMM "
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو SMP_ضربی_عدد 256
#ورنہ
#تعریف_کرو SMP_ضربی_عدد 128
#اگر_ختم

#اگر_تعریف_شدہ DYNAMIC_ARCH
خارجی حرف* gotoblas_corename(خالی);                                
#اگر_ختم                                                               

ساکن صحیح (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_غیر_تعریف_شدہ TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#ورنہ
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#اگر_ختم
};

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *SIDE, حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blas_arg_t دلائل;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم

  blasint info;
  صحیح side;
  صحیح uplo;
  صحیح unit;
  صحیح trans;
  صحیح nrowa;

  PRINT_DEBUG_NAME;

  دلائل.m = *M;
  دلائل.n = *N;

  دلائل.a = (خالی *)a;
  دلائل.b = (خالی *)b;

  دلائل.lda = *ldA;
  دلائل.ldb = *ldB;

  دلائل.beta = (خالی *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  اگر (side_arg  == 'L') side  = 0;
  اگر (side_arg  == 'R') side  = 1;

  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'R') trans = 2;
  اگر (trans_arg == 'C') trans = 3;

  اگر (diag_arg  == 'U') unit  = 0;
  اگر (diag_arg  == 'N') unit  = 1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  nrowa = دلائل.m;
  اگر (side & 1) nrowa = دلائل.n;

  info = 0;

  اگر (دلائل.ldb < MAX(1,دلائل.m)) info = 11;
  اگر (دلائل.lda < MAX(1,nrowa))  info =  9;
  اگر (دلائل.n < 0)               info =  6;
  اگر (دلائل.m < 0)               info =  5;
  اگر (unit < 0)                 info =  4;
  اگر (trans < 0)                info =  3;
  اگر (uplo  < 0)                info =  2;
  اگر (side  < 0)                info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME)-1);
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب,
	   شمار CBLAS_SIDE Side,  شمار CBLAS_UPLO Uplo,
	   شمار CBLAS_TRANSPOSE Trans, شمار CBLAS_DIAG Diag,
	   blasint m, blasint n,
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#ورنہ
	   خالی *valpha,
	   خالی *va, blasint lda,
	   خالی *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#اگر_ختم

  blas_arg_t دلائل;
  صحیح side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *بفر;
  XFLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم

  PRINT_DEBUG_CNAME;

  دلائل.a = (خالی *)a;
  دلائل.b = (خالی *)b;

  دلائل.lda = lda;
  دلائل.ldb = ldb;

#اگر_غیر_تعریف_شدہ COMPLEX
  دلائل.beta = (خالی *)&alpha;
#ورنہ
  دلائل.beta = (خالی *)alpha;
#اگر_ختم

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    دلائل.m = m;
    دلائل.n = n;

    اگر (Side == Cblas_بایاں)         side  = 0;
    اگر (Side == Cblas_دایاں)        side  = 1;

    اگر (Uplo == Cblas_بالائی)        uplo  = 0;
    اگر (Uplo == Cblas_زیریں)        uplo  = 1;

    اگر (Trans == CblasNoTrans)     trans = 0;
    اگر (Trans == CblasTrans)       trans = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (Trans == CblasConjNoTrans) trans = 0;
    اگر (Trans == CblasConjTrans)   trans = 1;
#ورنہ
    اگر (Trans == CblasConjNoTrans) trans = 2;
    اگر (Trans == CblasConjTrans)   trans = 3;
#اگر_ختم

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = دلائل.m;
    اگر (side & 1) nrowa = دلائل.n;

    اگر (دلائل.ldb < MAX(1,دلائل.m)) info = 11;
    اگر (دلائل.lda < MAX(1,nrowa))  info =  9;
    اگر (دلائل.n < 0)               info =  6;
    اگر (دلائل.m < 0)               info =  5;
    اگر (unit < 0)                 info =  4;
    اگر (trans < 0)                info =  3;
    اگر (uplo  < 0)                info =  2;
    اگر (side  < 0)                info =  1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    دلائل.m = n;
    دلائل.n = m;

    اگر (Side == Cblas_بایاں)         side  = 1;
    اگر (Side == Cblas_دایاں)        side  = 0;

    اگر (Uplo == Cblas_بالائی)        uplo  = 1;
    اگر (Uplo == Cblas_زیریں)        uplo  = 0;

    اگر (Trans == CblasNoTrans)     trans = 0;
    اگر (Trans == CblasTrans)       trans = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (Trans == CblasConjNoTrans) trans = 0;
    اگر (Trans == CblasConjTrans)   trans = 1;
#ورنہ
    اگر (Trans == CblasConjNoTrans) trans = 2;
    اگر (Trans == CblasConjTrans)   trans = 3;
#اگر_ختم

    اگر (Diag == CblasUnit)         unit  = 0;
    اگر (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = دلائل.m;
    اگر (side & 1) nrowa = دلائل.n;

    اگر (دلائل.ldb < MAX(1,دلائل.m)) info = 11;
    اگر (دلائل.lda < MAX(1,nrowa))  info =  9;
    اگر (دلائل.n < 0)               info =  6;
    اگر (دلائل.m < 0)               info =  5;
    اگر (unit < 0)                 info =  4;
    اگر (trans < 0)                info =  3;
    اگر (uplo  < 0)                info =  2;
    اگر (side  < 0)                info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

  اگر (دلائل.m == 0 || دلائل.n == 0) واپس;

#اگر defined(TRMM)
#اگر !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#اگر_ختم
)
#اگر_ختم
  اگر (ترتیب == Cblas_صف_Major && Diag == CblasNonUnit && Side == Cblas_بایاں && m == lda && n == ldb) {
    اگر (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_بالائی ? STRMM_براہ_راست_LNUN : STRMM_براہ_راست_LNLN)(m, n, alpha, a, lda, b, ldb);
      واپس;
    } ورنہ اگر (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_بالائی ? STRMM_براہ_راست_LTUN : STRMM_براہ_راست_LTLN)(m, n, alpha, a, lda, b, ldb);
      واپس;
    }
  }
#اگر_ختم
#اگر_ختم
#اگر_ختم



//end of the ifndef CBLAS ... else ... section	
#اگر_ختم

  اگر ((دلائل.m == 0) || (دلائل.n == 0)) واپس;
  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);

#اگر_تعریف_شدہ SMP
  منوال |= (trans << BLAS_TRANSA_سرکاؤ);
  منوال |= (side  << BLAS_RSIDE_سرکاؤ);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  اگر ( دلائل.m * دلائل.n < SMP_ضربی_عدد * GEMM_MULTITHREAD_حد)
	دلائل.nthreads = 1;
  ورنہ
	دلائل.nthreads = عدد_cpu_avail(3);
		

  اگر (دلائل.nthreads == 1) {
#اگر_ختم

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP
  } ورنہ {
    اگر (!side) {
      gemm_thread_n(منوال, &دلائل, تہی, تہی, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلائل.nthreads);
    } ورنہ {
      gemm_thread_m(منوال, &دلائل, تہی, تہی, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, دلائل.nthreads);
    }
  }
#اگر_ختم

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE,
		       (!side) ? دلائل.m * (دلائل.m + دلائل.n) : دلائل.n * (دلائل.m + دلائل.n),
		       (!side) ? دلائل.m * دلائل.m * دلائل.n : دلائل.m * دلائل.n * دلائل.n);

  IDEBUG_انتہا;

  واپس;
}

