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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_अपरिभाषितम् TRMM
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QTRSM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DTRSM "
#अन्यथा
#परिभाषय ERROR_NAME "STRSM "
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XTRSM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZTRSM "
#अन्यथा
#परिभाषय ERROR_NAME "CTRSM "
#यदि_अन्तः
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QTRMM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DTRMM "
#अन्यथा
#परिभाषय ERROR_NAME "STRMM "
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XTRMM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZTRMM "
#अन्यथा
#परिभाषय ERROR_NAME "CTRMM "
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि_अपरिभाषितम् COMPLEX
#परिभाषय SMP_गुणकः 256
#अन्यथा
#परिभाषय SMP_गुणकः 128
#यदि_अन्तः

#यदि_परिभाषितम् DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);                                
#यदि_अन्तः                                                               

स्थिर पूर्णाङ्क (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#यदि_अपरिभाषितम् TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#अन्यथा
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#यदि_अन्तः
};

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *SIDE, वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  वर्ण side_arg  = *SIDE;
  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blas_arg_t तर्काः;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  blasint info;
  पूर्णाङ्क side;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans;
  पूर्णाङ्क nrowa;

  PRINT_DEBUG_NAME;

  तर्काः.m = *M;
  तर्काः.n = *N;

  तर्काः.a = (शून्य *)a;
  तर्काः.b = (शून्य *)b;

  तर्काः.lda = *ldA;
  तर्काः.ldb = *ldB;

  तर्काः.beta = (शून्य *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  यदि (side_arg  == 'L') side  = 0;
  यदि (side_arg  == 'R') side  = 1;

  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'R') trans = 2;
  यदि (trans_arg == 'C') trans = 3;

  यदि (diag_arg  == 'U') unit  = 0;
  यदि (diag_arg  == 'N') unit  = 1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  nrowa = तर्काः.m;
  यदि (side & 1) nrowa = तर्काः.n;

  info = 0;

  यदि (तर्काः.ldb < MAX(1,तर्काः.m)) info = 11;
  यदि (तर्काः.lda < MAX(1,nrowa))  info =  9;
  यदि (तर्काः.n < 0)               info =  6;
  यदि (तर्काः.m < 0)               info =  5;
  यदि (unit < 0)                 info =  4;
  यदि (trans < 0)                info =  3;
  यदि (uplo  < 0)                info =  2;
  यदि (side  < 0)                info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME)-1);
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः,
	   गणना CBLAS_SIDE Side,  गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE Trans, गणना CBLAS_DIAG Diag,
	   blasint m, blasint n,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#यदि_अन्तः

  blas_arg_t तर्काः;
  पूर्णाङ्क side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *बफर;
  XFLOAT *sa, *sb;

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  पूर्णाङ्क बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  तर्काः.a = (शून्य *)a;
  तर्काः.b = (शून्य *)b;

  तर्काः.lda = lda;
  तर्काः.ldb = ldb;

#यदि_अपरिभाषितम् COMPLEX
  तर्काः.beta = (शून्य *)&alpha;
#अन्यथा
  तर्काः.beta = (शून्य *)alpha;
#यदि_अन्तः

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  यदि (क्रमः == CblasColMajor) {
    तर्काः.m = m;
    तर्काः.n = n;

    यदि (Side == Cblas_वामः)         side  = 0;
    यदि (Side == Cblas_दक्षिणः)        side  = 1;

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)        uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्)        uplo  = 1;

    यदि (Trans == CblasNoTrans)     trans = 0;
    यदि (Trans == CblasTrans)       trans = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    यदि (Trans == CblasConjNoTrans) trans = 2;
    यदि (Trans == CblasConjTrans)   trans = 3;
#यदि_अन्तः

    यदि (Diag == CblasUnit)          unit  = 0;
    यदि (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = तर्काः.m;
    यदि (side & 1) nrowa = तर्काः.n;

    यदि (तर्काः.ldb < MAX(1,तर्काः.m)) info = 11;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  9;
    यदि (तर्काः.n < 0)               info =  6;
    यदि (तर्काः.m < 0)               info =  5;
    यदि (unit < 0)                 info =  4;
    यदि (trans < 0)                info =  3;
    यदि (uplo  < 0)                info =  2;
    यदि (side  < 0)                info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    तर्काः.m = n;
    तर्काः.n = m;

    यदि (Side == Cblas_वामः)         side  = 1;
    यदि (Side == Cblas_दक्षिणः)        side  = 0;

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)        uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्)        uplo  = 0;

    यदि (Trans == CblasNoTrans)     trans = 0;
    यदि (Trans == CblasTrans)       trans = 1;
#यदि_अपरिभाषितम् COMPLEX
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    यदि (Trans == CblasConjNoTrans) trans = 2;
    यदि (Trans == CblasConjTrans)   trans = 3;
#यदि_अन्तः

    यदि (Diag == CblasUnit)         unit  = 0;
    यदि (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = तर्काः.m;
    यदि (side & 1) nrowa = तर्काः.n;

    यदि (तर्काः.ldb < MAX(1,तर्काः.m)) info = 11;
    यदि (तर्काः.lda < MAX(1,nrowa))  info =  9;
    यदि (तर्काः.n < 0)               info =  6;
    यदि (तर्काः.m < 0)               info =  5;
    यदि (unit < 0)                 info =  4;
    यदि (trans < 0)                info =  3;
    यदि (uplo  < 0)                info =  2;
    यदि (side  < 0)                info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (तर्काः.m == 0 || तर्काः.n == 0) प्रत्यागम;

#यदि defined(TRMM)
#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः
  यदि (क्रमः == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_वामः && m == lda && n == ldb) {
    यदि (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_ऊर्ध्वस्थम् ? STRMM_प्रत्यक्षम्_LNUN : STRMM_प्रत्यक्षम्_LNLN)(m, n, alpha, a, lda, b, ldb);
      प्रत्यागम;
    } अन्यथा यदि (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_ऊर्ध्वस्थम् ? STRMM_प्रत्यक्षम्_LTUN : STRMM_प्रत्यक्षम्_LTLN)(m, n, alpha, a, lda, b, ldb);
      प्रत्यागम;
    }
  }
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः



//end of the ifndef CBLAS ... else ... section	
#यदि_अन्तः

  यदि ((तर्काः.m == 0) || (तर्काः.n == 0)) प्रत्यागम;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  बफर = (FLOAT *)blas_स्मृतिः_alloc(0);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);

#यदि_परिभाषितम् SMP
  बहुलकम् |= (trans << BLAS_TRANSA_सरणम्);
  बहुलकम् |= (side  << BLAS_RSIDE_सरणम्);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  यदि ( तर्काः.m * तर्काः.n < SMP_गुणकः * GEMM_MULTITHREAD_देहली)
	तर्काः.nthreads = 1;
  अन्यथा
	तर्काः.nthreads = अङ्कः_cpu_avail(3);
		

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP
  } अन्यथा {
    यदि (!side) {
      gemm_thread_n(बहुलकम्, &तर्काः, रिक्त, रिक्त, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, तर्काः.nthreads);
    } अन्यथा {
      gemm_thread_m(बहुलकम्, &तर्काः, रिक्त, रिक्त, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, तर्काः.nthreads);
    }
  }
#यदि_अन्तः

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE,
		       (!side) ? तर्काः.m * (तर्काः.m + तर्काः.n) : तर्काः.n * (तर्काः.m + तर्काः.n),
		       (!side) ? तर्काः.m * तर्काः.m * तर्काः.n : तर्काः.m * तर्काः.n * तर्काः.n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}

