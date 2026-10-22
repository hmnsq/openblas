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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ TRMM
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QTRSM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DTRSM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "STRSM "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XTRSM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZTRSM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CTRSM "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QTRMM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DTRMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "STRMM "
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XTRMM "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZTRMM "
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CTRMM "
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಗುಣಕ 256
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMP_ಗುಣಕ 128
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DYNAMIC_ARCH
ಬಾಹ್ಯ ಅಕ್ಷರ* gotoblas_corename(ಶೂನ್ಯ);                                
#ಆದರೆ_ಅಂತ್ಯ                                                               

ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#ಇಲ್ಲದಿದ್ದರೆ
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#ಆದರೆ_ಅಂತ್ಯ
};

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(ಅಕ್ಷರ *SIDE, ಅಕ್ಷರ *UPLO, ಅಕ್ಷರ *TRANS, ಅಕ್ಷರ *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  ಅಕ್ಷರ side_arg  = *SIDE;
  ಅಕ್ಷರ uplo_arg  = *UPLO;
  ಅಕ್ಷರ trans_arg = *TRANS;
  ಅಕ್ಷರ diag_arg  = *DIAG;

  blas_arg_t ವಾದಗಳು;

  FLOAT *ಬಫರ್;
  FLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

  blasint info;
  ಪೂರ್ಣಾಂಕ side;
  ಪೂರ್ಣಾಂಕ uplo;
  ಪೂರ್ಣಾಂಕ unit;
  ಪೂರ್ಣಾಂಕ trans;
  ಪೂರ್ಣಾಂಕ nrowa;

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.m = *M;
  ವಾದಗಳು.n = *N;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;

  ವಾದಗಳು.lda = *ldA;
  ವಾದಗಳು.ldb = *ldB;

  ವಾದಗಳು.beta = (ಶೂನ್ಯ *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  ಆದರೆ (side_arg  == 'L') side  = 0;
  ಆದರೆ (side_arg  == 'R') side  = 1;

  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'T') trans = 1;
  ಆದರೆ (trans_arg == 'R') trans = 2;
  ಆದರೆ (trans_arg == 'C') trans = 3;

  ಆದರೆ (diag_arg  == 'U') unit  = 0;
  ಆದರೆ (diag_arg  == 'N') unit  = 1;

  ಆದರೆ (uplo_arg  == 'U') uplo  = 0;
  ಆದರೆ (uplo_arg  == 'L') uplo  = 1;

  nrowa = ವಾದಗಳು.m;
  ಆದರೆ (side & 1) nrowa = ವಾದಗಳು.n;

  info = 0;

  ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,ವಾದಗಳು.m)) info = 11;
  ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  9;
  ಆದರೆ (ವಾದಗಳು.n < 0)               info =  6;
  ಆದರೆ (ವಾದಗಳು.m < 0)               info =  5;
  ಆದರೆ (unit < 0)                 info =  4;
  ಆದರೆ (trans < 0)                info =  3;
  ಆದರೆ (uplo  < 0)                info =  2;
  ಆದರೆ (side  < 0)                info =  1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME)-1);
    ಹಿಂತಿರುಗಿ;
  }

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ,
	   ಎಣಿಕೆ CBLAS_SIDE Side,  ಎಣಿಕೆ CBLAS_UPLO Uplo,
	   ಎಣಿಕೆ CBLAS_TRANSPOSE Trans, ಎಣಿಕೆ CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ *valpha,
	   ಶೂನ್ಯ *va, blasint lda,
	   ಶೂನ್ಯ *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#ಆದರೆ_ಅಂತ್ಯ

  blas_arg_t ವಾದಗಳು;
  ಪೂರ್ಣಾಂಕ side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *ಬಫರ್;
  XFLOAT *sa, *sb;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;

  ವಾದಗಳು.a = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.b = (ಶೂನ್ಯ *)b;

  ವಾದಗಳು.lda = lda;
  ವಾದಗಳು.ldb = ldb;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ವಾದಗಳು.beta = (ಶೂನ್ಯ *)&alpha;
#ಇಲ್ಲದಿದ್ದರೆ
  ವಾದಗಳು.beta = (ಶೂನ್ಯ *)alpha;
#ಆದರೆ_ಅಂತ್ಯ

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
    ವಾದಗಳು.m = m;
    ವಾದಗಳು.n = n;

    ಆದರೆ (Side == Cblas_ಎಡ)         side  = 0;
    ಆದರೆ (Side == Cblas_ಬಲ)        side  = 1;

    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)        uplo  = 0;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)        uplo  = 1;

    ಆದರೆ (Trans == CblasNoTrans)     trans = 0;
    ಆದರೆ (Trans == CblasTrans)       trans = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 0;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 2;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 3;
#ಆದರೆ_ಅಂತ್ಯ

    ಆದರೆ (Diag == CblasUnit)          unit  = 0;
    ಆದರೆ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = ವಾದಗಳು.m;
    ಆದರೆ (side & 1) nrowa = ವಾದಗಳು.n;

    ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,ವಾದಗಳು.m)) info = 11;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  9;
    ಆದರೆ (ವಾದಗಳು.n < 0)               info =  6;
    ಆದರೆ (ವಾದಗಳು.m < 0)               info =  5;
    ಆದರೆ (unit < 0)                 info =  4;
    ಆದರೆ (trans < 0)                info =  3;
    ಆದರೆ (uplo  < 0)                info =  2;
    ಆದರೆ (side  < 0)                info =  1;
  }

  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
    ವಾದಗಳು.m = n;
    ವಾದಗಳು.n = m;

    ಆದರೆ (Side == Cblas_ಎಡ)         side  = 1;
    ಆದರೆ (Side == Cblas_ಬಲ)        side  = 0;

    ಆದರೆ (Uplo == Cblas_ಊರ್ಧ್ವ)        uplo  = 1;
    ಆದರೆ (Uplo == Cblas_ಅಧೋ)        uplo  = 0;

    ಆದರೆ (Trans == CblasNoTrans)     trans = 0;
    ಆದರೆ (Trans == CblasTrans)       trans = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 0;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 1;
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (Trans == CblasConjNoTrans) trans = 2;
    ಆದರೆ (Trans == CblasConjTrans)   trans = 3;
#ಆದರೆ_ಅಂತ್ಯ

    ಆದರೆ (Diag == CblasUnit)         unit  = 0;
    ಆದರೆ (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = ವಾದಗಳು.m;
    ಆದರೆ (side & 1) nrowa = ವಾದಗಳು.n;

    ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,ವಾದಗಳು.m)) info = 11;
    ಆದರೆ (ವಾದಗಳು.lda < MAX(1,nrowa))  info =  9;
    ಆದರೆ (ವಾದಗಳು.n < 0)               info =  6;
    ಆದರೆ (ವಾದಗಳು.m < 0)               info =  5;
    ಆದರೆ (unit < 0)                 info =  4;
    ಆದರೆ (trans < 0)                info =  3;
    ಆದರೆ (uplo  < 0)                info =  2;
    ಆದರೆ (side  < 0)                info =  1;
  }

  ಆದರೆ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

  ಆದರೆ (ವಾದಗಳು.m == 0 || ವಾದಗಳು.n == 0) ಹಿಂತಿರುಗಿ;

#ಆದರೆ defined(TRMM)
#ಆದರೆ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಆದರೆ defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_ನೇರ)||defined(DYNAMIC_ARCH))
#ಆದರೆ defined(DYNAMIC_ARCH)
ಆದರೆ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಆದರೆ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ಆದರೆ_ಅಂತ್ಯ
)
#ಆದರೆ_ಅಂತ್ಯ
  ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major && Diag == CblasNonUnit && Side == Cblas_ಎಡ && m == lda && n == ldb) {
    ಆದರೆ (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_ಊರ್ಧ್ವ ? STRMM_ನೇರ_LNUN : STRMM_ನೇರ_LNLN)(m, n, alpha, a, lda, b, ldb);
      ಹಿಂತಿರುಗಿ;
    } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_ಊರ್ಧ್ವ ? STRMM_ನೇರ_LTUN : STRMM_ನೇರ_LTLN)(m, n, alpha, a, lda, b, ldb);
      ಹಿಂತಿರುಗಿ;
    }
  }
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ



//end of the ifndef CBLAS ... else ... section	
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ ((ವಾದಗಳು.m == 0) || (ವಾದಗಳು.n == 0)) ಹಿಂತಿರುಗಿ;
  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(0);
  ಆದರೆ (!ಬಫರ್) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
    ಹಿಂತಿರುಗಿ;
  }

  sa = (FLOAT *)((BLASLONG)ಬಫರ್ + GEMM_ಸ್ಥಾನಾಂತರ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಬಹುಲಕ |= (trans << BLAS_TRANSA_ಸರಿಸು);
  ಬಹುಲಕ |= (side  << BLAS_RSIDE_ಸರಿಸು);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  ಆದರೆ ( ವಾದಗಳು.m * ವಾದಗಳು.n < SMP_ಗುಣಕ * GEMM_MULTITHREAD_ಮಿತಿ)
	ವಾದಗಳು.nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ
	ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(3);
		

  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {
    ಆದರೆ (!side) {
      gemm_thread_n(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, ವಾದಗಳು.nthreads);
    } ಇಲ್ಲದಿದ್ದರೆ {
      gemm_thread_m(ಬಹುಲಕ, &ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, ವಾದಗಳು.nthreads);
    }
  }
#ಆದರೆ_ಅಂತ್ಯ

  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE,
		       (!side) ? ವಾದಗಳು.m * (ವಾದಗಳು.m + ವಾದಗಳು.n) : ವಾದಗಳು.n * (ವಾದಗಳು.m + ವಾದಗಳು.n),
		       (!side) ? ವಾದಗಳು.m * ವಾದಗಳು.m * ವಾದಗಳು.n : ವಾದಗಳು.m * ವಾದಗಳು.n * ವಾದಗಳು.n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}

