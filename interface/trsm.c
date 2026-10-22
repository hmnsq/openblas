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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <ctype.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ TRMM
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QTRSM "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DTRSM "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "STRSM "
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTRSM "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTRSM "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTRSM "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QTRMM "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DTRMM "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "STRMM "
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTRMM "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTRMM "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTRMM "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਗੁਣਕ 256
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਗੁਣਕ 128
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
ਬਾਹਰੀ ਅੱਖਰ* gotoblas_corename(ਖਾਲੀ);                                
#ਜੇ_ਅੰਤ                                                               

ਸਥਿਰ ਪੂਰਨ (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#ਨਹੀਂਤਾਂ
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#ਜੇ_ਅੰਤ
};

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *SIDE, ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANS, ਅੱਖਰ *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  ਅੱਖਰ side_arg  = *SIDE;
  ਅੱਖਰ uplo_arg  = *UPLO;
  ਅੱਖਰ trans_arg = *TRANS;
  ਅੱਖਰ diag_arg  = *DIAG;

  blas_arg_t ਦਲੀਲਾਂ;

  FLOAT *ਬਫਰ;
  FLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  blasint info;
  ਪੂਰਨ side;
  ਪੂਰਨ uplo;
  ਪੂਰਨ unit;
  ਪੂਰਨ trans;
  ਪੂਰਨ nrowa;

  PRINT_DEBUG_NAME;

  ਦਲੀਲਾਂ.m = *M;
  ਦਲੀਲਾਂ.n = *N;

  ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)b;

  ਦਲੀਲਾਂ.lda = *ldA;
  ਦਲੀਲਾਂ.ldb = *ldB;

  ਦਲੀਲਾਂ.beta = (ਖਾਲੀ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  ਜੇ (side_arg  == 'L') side  = 0;
  ਜੇ (side_arg  == 'R') side  = 1;

  ਜੇ (trans_arg == 'N') trans = 0;
  ਜੇ (trans_arg == 'T') trans = 1;
  ਜੇ (trans_arg == 'R') trans = 2;
  ਜੇ (trans_arg == 'C') trans = 3;

  ਜੇ (diag_arg  == 'U') unit  = 0;
  ਜੇ (diag_arg  == 'N') unit  = 1;

  ਜੇ (uplo_arg  == 'U') uplo  = 0;
  ਜੇ (uplo_arg  == 'L') uplo  = 1;

  nrowa = ਦਲੀਲਾਂ.m;
  ਜੇ (side & 1) nrowa = ਦਲੀਲਾਂ.n;

  info = 0;

  ਜੇ (ਦਲੀਲਾਂ.ldb < MAX(1,ਦਲੀਲਾਂ.m)) info = 11;
  ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  9;
  ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  6;
  ਜੇ (ਦਲੀਲਾਂ.m < 0)               info =  5;
  ਜੇ (unit < 0)                 info =  4;
  ਜੇ (trans < 0)                info =  3;
  ਜੇ (uplo  < 0)                info =  2;
  ਜੇ (side  < 0)                info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME)-1);
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ,
	   ਗਿਣਤੀ CBLAS_SIDE Side,  ਗਿਣਤੀ CBLAS_UPLO Uplo,
	   ਗਿਣਤੀ CBLAS_TRANSPOSE Trans, ਗਿਣਤੀ CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *valpha,
	   ਖਾਲੀ *va, blasint lda,
	   ਖਾਲੀ *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#ਜੇ_ਅੰਤ

  blas_arg_t ਦਲੀਲਾਂ;
  ਪੂਰਨ side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *ਬਫਰ;
  XFLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.b = (ਖਾਲੀ *)b;

  ਦਲੀਲਾਂ.lda = lda;
  ਦਲੀਲਾਂ.ldb = ldb;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਦਲੀਲਾਂ.beta = (ਖਾਲੀ *)&alpha;
#ਨਹੀਂਤਾਂ
  ਦਲੀਲਾਂ.beta = (ਖਾਲੀ *)alpha;
#ਜੇ_ਅੰਤ

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਦਲੀਲਾਂ.m = m;
    ਦਲੀਲਾਂ.n = n;

    ਜੇ (Side == Cblas_ਖੱਬਾ)         side  = 0;
    ਜੇ (Side == Cblas_ਸੱਜਾ)        side  = 1;

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)        uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)        uplo  = 1;

    ਜੇ (Trans == CblasNoTrans)     trans = 0;
    ਜੇ (Trans == CblasTrans)       trans = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (Trans == CblasConjNoTrans) trans = 0;
    ਜੇ (Trans == CblasConjTrans)   trans = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (Trans == CblasConjNoTrans) trans = 2;
    ਜੇ (Trans == CblasConjTrans)   trans = 3;
#ਜੇ_ਅੰਤ

    ਜੇ (Diag == CblasUnit)          unit  = 0;
    ਜੇ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = ਦਲੀਲਾਂ.m;
    ਜੇ (side & 1) nrowa = ਦਲੀਲਾਂ.n;

    ਜੇ (ਦਲੀਲਾਂ.ldb < MAX(1,ਦਲੀਲਾਂ.m)) info = 11;
    ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  9;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  6;
    ਜੇ (ਦਲੀਲਾਂ.m < 0)               info =  5;
    ਜੇ (unit < 0)                 info =  4;
    ਜੇ (trans < 0)                info =  3;
    ਜੇ (uplo  < 0)                info =  2;
    ਜੇ (side  < 0)                info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਦਲੀਲਾਂ.m = n;
    ਦਲੀਲਾਂ.n = m;

    ਜੇ (Side == Cblas_ਖੱਬਾ)         side  = 1;
    ਜੇ (Side == Cblas_ਸੱਜਾ)        side  = 0;

    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ)        uplo  = 1;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ)        uplo  = 0;

    ਜੇ (Trans == CblasNoTrans)     trans = 0;
    ਜੇ (Trans == CblasTrans)       trans = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (Trans == CblasConjNoTrans) trans = 0;
    ਜੇ (Trans == CblasConjTrans)   trans = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (Trans == CblasConjNoTrans) trans = 2;
    ਜੇ (Trans == CblasConjTrans)   trans = 3;
#ਜੇ_ਅੰਤ

    ਜੇ (Diag == CblasUnit)         unit  = 0;
    ਜੇ (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = ਦਲੀਲਾਂ.m;
    ਜੇ (side & 1) nrowa = ਦਲੀਲਾਂ.n;

    ਜੇ (ਦਲੀਲਾਂ.ldb < MAX(1,ਦਲੀਲਾਂ.m)) info = 11;
    ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  9;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  6;
    ਜੇ (ਦਲੀਲਾਂ.m < 0)               info =  5;
    ਜੇ (unit < 0)                 info =  4;
    ਜੇ (trans < 0)                info =  3;
    ਜੇ (uplo  < 0)                info =  2;
    ਜੇ (side  < 0)                info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

  ਜੇ (ਦਲੀਲਾਂ.m == 0 || ਦਲੀਲਾਂ.n == 0) ਵਾਪਸ;

#ਜੇ defined(TRMM)
#ਜੇ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ਜੇ defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_ਸਿੱਧਾ)||defined(DYNAMIC_ARCH))
#ਜੇ defined(DYNAMIC_ARCH)
ਜੇ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ਜੇ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ਜੇ_ਅੰਤ
)
#ਜੇ_ਅੰਤ
  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major && Diag == CblasNonUnit && Side == Cblas_ਖੱਬਾ && m == lda && n == ldb) {
    ਜੇ (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_ਉੱਤਲਾ ? STRMM_ਸਿੱਧਾ_LNUN : STRMM_ਸਿੱਧਾ_LNLN)(m, n, alpha, a, lda, b, ldb);
      ਵਾਪਸ;
    } ਨਹੀਂਤਾਂ ਜੇ (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_ਉੱਤਲਾ ? STRMM_ਸਿੱਧਾ_LTUN : STRMM_ਸਿੱਧਾ_LTLN)(m, n, alpha, a, lda, b, ldb);
      ਵਾਪਸ;
    }
  }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ



//end of the ifndef CBLAS ... else ... section	
#ਜੇ_ਅੰਤ

  ਜੇ ((ਦਲੀਲਾਂ.m == 0) || (ਦਲੀਲਾਂ.n == 0)) ਵਾਪਸ;
  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(0);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

  sa = (FLOAT *)((BLASLONG)ਬਫਰ + GEMM_ਵਿੱਥ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਬਹੁਲਕ |= (trans << BLAS_TRANSA_ਖਿਸਕਾਓ);
  ਬਹੁਲਕ |= (side  << BLAS_RSIDE_ਖਿਸਕਾਓ);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  ਜੇ ( ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n < SMP_ਗੁਣਕ * GEMM_MULTITHREAD_ਹੱਦ)
	ਦਲੀਲਾਂ.nthreads = 1;
  ਨਹੀਂਤਾਂ
	ਦਲੀਲਾਂ.nthreads = ਅੰਕ_cpu_avail(3);
		

  ਜੇ (ਦਲੀਲਾਂ.nthreads == 1) {
#ਜੇ_ਅੰਤ

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {
    ਜੇ (!side) {
      gemm_thread_n(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, ਦਲੀਲਾਂ.nthreads);
    } ਨਹੀਂਤਾਂ {
      gemm_thread_m(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, ਦਲੀਲਾਂ.nthreads);
    }
  }
#ਜੇ_ਅੰਤ

  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE * COMPSIZE,
		       (!side) ? ਦਲੀਲਾਂ.m * (ਦਲੀਲਾਂ.m + ਦਲੀਲਾਂ.n) : ਦਲੀਲਾਂ.n * (ਦਲੀਲਾਂ.m + ਦਲੀਲਾਂ.n),
		       (!side) ? ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n : ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}

