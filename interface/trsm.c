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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <ctype.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే TRMM
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QTRSM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DTRSM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "STRSM "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XTRSM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZTRSM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CTRSM "
#అయితే_ముగింపు
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QTRMM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DTRMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "STRMM "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XTRMM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZTRMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CTRMM "
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు SMP_గుణకం 256
#లేకపోతే
#నిర్వచించు SMP_గుణకం 128
#అయితే_ముగింపు

#నిర్వచితమైతే DYNAMIC_ARCH
బాహ్య అక్షరం* gotoblas_corename(శూన్య);                                
#అయితే_ముగింపు                                                               

స్థిర పూర్ణాంకం (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#నిర్వచితం_కాకపోతే TRMM
  TRSM_LNUU, TRSM_LNUN, TRSM_LNLU, TRSM_LNLN,
  TRSM_LTUU, TRSM_LTUN, TRSM_LTLU, TRSM_LTLN,
  TRSM_LRUU, TRSM_LRUN, TRSM_LRLU, TRSM_LRLN,
  TRSM_LCUU, TRSM_LCUN, TRSM_LCLU, TRSM_LCLN,
  TRSM_RNUU, TRSM_RNUN, TRSM_RNLU, TRSM_RNLN,
  TRSM_RTUU, TRSM_RTUN, TRSM_RTLU, TRSM_RTLN,
  TRSM_RRUU, TRSM_RRUN, TRSM_RRLU, TRSM_RRLN,
  TRSM_RCUU, TRSM_RCUN, TRSM_RCLU, TRSM_RCLN,
#లేకపోతే
  TRMM_LNUU, TRMM_LNUN, TRMM_LNLU, TRMM_LNLN,
  TRMM_LTUU, TRMM_LTUN, TRMM_LTLU, TRMM_LTLN,
  TRMM_LRUU, TRMM_LRUN, TRMM_LRLU, TRMM_LRLN,
  TRMM_LCUU, TRMM_LCUN, TRMM_LCLU, TRMM_LCLN,
  TRMM_RNUU, TRMM_RNUN, TRMM_RNLU, TRMM_RNLN,
  TRMM_RTUU, TRMM_RTUN, TRMM_RTLU, TRMM_RTLN,
  TRMM_RRUU, TRMM_RRUN, TRMM_RRLU, TRMM_RRLN,
  TRMM_RCUU, TRMM_RCUN, TRMM_RCLU, TRMM_RCLN,
#అయితే_ముగింపు
};

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *SIDE, అక్షరం *UPLO, అక్షరం *TRANS, అక్షరం *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  అక్షరం side_arg  = *SIDE;
  అక్షరం uplo_arg  = *UPLO;
  అక్షరం trans_arg = *TRANS;
  అక్షరం diag_arg  = *DIAG;

  blas_arg_t వాదనలు;

  FLOAT *బఫరు;
  FLOAT *sa, *sb;

#నిర్వచితమైతే SMP
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

  blasint info;
  పూర్ణాంకం side;
  పూర్ణాంకం uplo;
  పూర్ణాంకం unit;
  పూర్ణాంకం trans;
  పూర్ణాంకం nrowa;

  PRINT_DEBUG_NAME;

  వాదనలు.m = *M;
  వాదనలు.n = *N;

  వాదనలు.a = (శూన్య *)a;
  వాదనలు.b = (శూన్య *)b;

  వాదనలు.lda = *ldA;
  వాదనలు.ldb = *ldB;

  వాదనలు.beta = (శూన్య *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  ఒకవేళ (side_arg  == 'L') side  = 0;
  ఒకవేళ (side_arg  == 'R') side  = 1;

  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
  ఒకవేళ (trans_arg == 'R') trans = 2;
  ఒకవేళ (trans_arg == 'C') trans = 3;

  ఒకవేళ (diag_arg  == 'U') unit  = 0;
  ఒకవేళ (diag_arg  == 'N') unit  = 1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;

  nrowa = వాదనలు.m;
  ఒకవేళ (side & 1) nrowa = వాదనలు.n;

  info = 0;

  ఒకవేళ (వాదనలు.ldb < MAX(1,వాదనలు.m)) info = 11;
  ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  9;
  ఒకవేళ (వాదనలు.n < 0)               info =  6;
  ఒకవేళ (వాదనలు.m < 0)               info =  5;
  ఒకవేళ (unit < 0)                 info =  4;
  ఒకవేళ (trans < 0)                info =  3;
  ఒకవేళ (uplo  < 0)                info =  2;
  ఒకవేళ (side  < 0)                info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME)-1);
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం,
	   గణన CBLAS_SIDE Side,  గణన CBLAS_UPLO Uplo,
	   గణన CBLAS_TRANSPOSE Trans, గణన CBLAS_DIAG Diag,
	   blasint m, blasint n,
#నిర్వచితం_కాకపోతే COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#లేకపోతే
	   శూన్య *valpha,
	   శూన్య *va, blasint lda,
	   శూన్య *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#అయితే_ముగింపు

  blas_arg_t వాదనలు;
  పూర్ణాంకం side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *బఫరు;
  XFLOAT *sa, *sb;

#నిర్వచితమైతే SMP
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  పూర్ణాంకం బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  పూర్ణాంకం బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  పూర్ణాంకం బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  వాదనలు.a = (శూన్య *)a;
  వాదనలు.b = (శూన్య *)b;

  వాదనలు.lda = lda;
  వాదనలు.ldb = ldb;

#నిర్వచితం_కాకపోతే COMPLEX
  వాదనలు.beta = (శూన్య *)&alpha;
#లేకపోతే
  వాదనలు.beta = (శూన్య *)alpha;
#అయితే_ముగింపు

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    వాదనలు.m = m;
    వాదనలు.n = n;

    ఒకవేళ (Side == Cblas_ఎడమ)         side  = 0;
    ఒకవేళ (Side == Cblas_కుడి)        side  = 1;

    ఒకవేళ (Uplo == Cblas_పై)        uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది)        uplo  = 1;

    ఒకవేళ (Trans == CblasNoTrans)     trans = 0;
    ఒకవేళ (Trans == CblasTrans)       trans = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 0;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 1;
#లేకపోతే
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 2;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 3;
#అయితే_ముగింపు

    ఒకవేళ (Diag == CblasUnit)          unit  = 0;
    ఒకవేళ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = వాదనలు.m;
    ఒకవేళ (side & 1) nrowa = వాదనలు.n;

    ఒకవేళ (వాదనలు.ldb < MAX(1,వాదనలు.m)) info = 11;
    ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  9;
    ఒకవేళ (వాదనలు.n < 0)               info =  6;
    ఒకవేళ (వాదనలు.m < 0)               info =  5;
    ఒకవేళ (unit < 0)                 info =  4;
    ఒకవేళ (trans < 0)                info =  3;
    ఒకవేళ (uplo  < 0)                info =  2;
    ఒకవేళ (side  < 0)                info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    వాదనలు.m = n;
    వాదనలు.n = m;

    ఒకవేళ (Side == Cblas_ఎడమ)         side  = 1;
    ఒకవేళ (Side == Cblas_కుడి)        side  = 0;

    ఒకవేళ (Uplo == Cblas_పై)        uplo  = 1;
    ఒకవేళ (Uplo == Cblas_క్రింది)        uplo  = 0;

    ఒకవేళ (Trans == CblasNoTrans)     trans = 0;
    ఒకవేళ (Trans == CblasTrans)       trans = 1;
#నిర్వచితం_కాకపోతే COMPLEX
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 0;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 1;
#లేకపోతే
    ఒకవేళ (Trans == CblasConjNoTrans) trans = 2;
    ఒకవేళ (Trans == CblasConjTrans)   trans = 3;
#అయితే_ముగింపు

    ఒకవేళ (Diag == CblasUnit)         unit  = 0;
    ఒకవేళ (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = వాదనలు.m;
    ఒకవేళ (side & 1) nrowa = వాదనలు.n;

    ఒకవేళ (వాదనలు.ldb < MAX(1,వాదనలు.m)) info = 11;
    ఒకవేళ (వాదనలు.lda < MAX(1,nrowa))  info =  9;
    ఒకవేళ (వాదనలు.n < 0)               info =  6;
    ఒకవేళ (వాదనలు.m < 0)               info =  5;
    ఒకవేళ (unit < 0)                 info =  4;
    ఒకవేళ (trans < 0)                info =  3;
    ఒకవేళ (uplo  < 0)                info =  2;
    ఒకవేళ (side  < 0)                info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

  ఒకవేళ (వాదనలు.m == 0 || వాదనలు.n == 0) తిరిగి;

#ఒకవేళ defined(TRMM)
#ఒకవేళ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ఒకవేళ defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
ఒకవేళ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ఒకవేళ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#అయితే_ముగింపు
)
#అయితే_ముగింపు
  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && Diag == CblasNonUnit && Side == Cblas_ఎడమ && m == lda && n == ldb) {
    ఒకవేళ (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == Cblas_పై ? STRMM_ప్రత్యక్ష_LNUN : STRMM_ప్రత్యక్ష_LNLN)(m, n, alpha, a, lda, b, ldb);
      తిరిగి;
    } లేకపోతే ఒకవేళ (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == Cblas_పై ? STRMM_ప్రత్యక్ష_LTUN : STRMM_ప్రత్యక్ష_LTLN)(m, n, alpha, a, lda, b, ldb);
      తిరిగి;
    }
  }
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు



//end of the ifndef CBLAS ... else ... section	
#అయితే_ముగింపు

  ఒకవేళ ((వాదనలు.m == 0) || (వాదనలు.n == 0)) తిరిగి;
  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(0);
  ఒకవేళ (!బఫరు) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

  sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);

#నిర్వచితమైతే SMP
  బహుళకం |= (trans << BLAS_TRANSA_జరుపు);
  బహుళకం |= (side  << BLAS_RSIDE_జరుపు);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  ఒకవేళ ( వాదనలు.m * వాదనలు.n < SMP_గుణకం * GEMM_MULTITHREAD_పరిమితి)
	వాదనలు.nthreads = 1;
  లేకపోతే
	వాదనలు.nthreads = అంకె_cpu_avail(3);
		

  ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP
  } లేకపోతే {
    ఒకవేళ (!side) {
      gemm_thread_n(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, వాదనలు.nthreads);
    } లేకపోతే {
      gemm_thread_m(బహుళకం, &వాదనలు, ఖాళీ, ఖాళీ, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, వాదనలు.nthreads);
    }
  }
#అయితే_ముగింపు

  blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE,
		       (!side) ? వాదనలు.m * (వాదనలు.m + వాదనలు.n) : వాదనలు.n * (వాదనలు.m + వాదనలు.n),
		       (!side) ? వాదనలు.m * వాదనలు.m * వాదనలు.n : వాదనలు.m * వాదనలు.n * వాదనలు.n);

  IDEBUG_చివర;

  తిరిగి;
}

