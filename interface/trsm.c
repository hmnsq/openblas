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

#समावेश <stdio.h>
#समावेश <ctype.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_अघोषित TRMM
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QTRSM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DTRSM "
#अन्यथा
#घोषणा ERROR_NAME "STRSM "
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XTRSM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZTRSM "
#अन्यथा
#घोषणा ERROR_NAME "CTRSM "
#पूर्ण_अगर
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QTRMM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DTRMM "
#अन्यथा
#घोषणा ERROR_NAME "STRMM "
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XTRMM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZTRMM "
#अन्यथा
#घोषणा ERROR_NAME "CTRMM "
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर_अघोषित COMPLEX
#घोषणा SMP_FACTOR 256
#अन्यथा
#घोषणा SMP_FACTOR 128
#पूर्ण_अगर

#अगर_घोषित DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(व्योम);                                
#पूर्ण_अगर                                                               

अटल पूर्णांक (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#अगर_अघोषित TRMM
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
#पूर्ण_अगर
};

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *SIDE, अक्षर *UPLO, अक्षर *TRANS, अक्षर *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  अक्षर side_arg  = *SIDE;
  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;
  अक्षर diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर_घोषित SMP
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

  blasint info;
  पूर्णांक side;
  पूर्णांक uplo;
  पूर्णांक unit;
  पूर्णांक trans;
  पूर्णांक nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (व्योम *)a;
  args.b = (व्योम *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (व्योम *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  अगर (side_arg  == 'L') side  = 0;
  अगर (side_arg  == 'R') side  = 1;

  अगर (trans_arg == 'N') trans = 0;
  अगर (trans_arg == 'T') trans = 1;
  अगर (trans_arg == 'R') trans = 2;
  अगर (trans_arg == 'C') trans = 3;

  अगर (diag_arg  == 'U') unit  = 0;
  अगर (diag_arg  == 'N') unit  = 1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  अगर (side & 1) nrowa = args.n;

  info = 0;

  अगर (args.ldb < MAX(1,args.m)) info = 11;
  अगर (args.lda < MAX(1,nrowa))  info =  9;
  अगर (args.n < 0)               info =  6;
  अगर (args.m < 0)               info =  5;
  अगर (unit < 0)                 info =  4;
  अगर (trans < 0)                info =  3;
  अगर (uplo  < 0)                info =  2;
  अगर (side  < 0)                info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME)-1);
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order,
	   क्रमागत CBLAS_SIDE Side,  क्रमागत CBLAS_UPLO Uplo,
	   क्रमागत CBLAS_TRANSPOSE Trans, क्रमागत CBLAS_DIAG Diag,
	   blasint m, blasint n,
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb) {
#अन्यथा
	   व्योम *valpha,
	   व्योम *va, blasint lda,
	   व्योम *vb, blasint ldb) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
#पूर्ण_अगर

  blas_arg_t args;
  पूर्णांक side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *बफर;
  XFLOAT *sa, *sb;

#अगर_घोषित SMP
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  args.a = (व्योम *)a;
  args.b = (व्योम *)b;

  args.lda = lda;
  args.ldb = ldb;

#अगर_अघोषित COMPLEX
  args.beta = (व्योम *)&alpha;
#अन्यथा
  args.beta = (व्योम *)alpha;
#पूर्ण_अगर

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  अगर (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    अगर (Side == Cblas_बायाँ)         side  = 0;
    अगर (Side == Cblas_दायाँ)        side  = 1;

    अगर (Uplo == CblasUpper)        uplo  = 0;
    अगर (Uplo == CblasLower)        uplo  = 1;

    अगर (Trans == CblasNoTrans)     trans = 0;
    अगर (Trans == CblasTrans)       trans = 1;
#अगर_अघोषित COMPLEX
    अगर (Trans == CblasConjNoTrans) trans = 0;
    अगर (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    अगर (Trans == CblasConjNoTrans) trans = 2;
    अगर (Trans == CblasConjTrans)   trans = 3;
#पूर्ण_अगर

    अगर (Diag == CblasUnit)          unit  = 0;
    अगर (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    अगर (side & 1) nrowa = args.n;

    अगर (args.ldb < MAX(1,args.m)) info = 11;
    अगर (args.lda < MAX(1,nrowa))  info =  9;
    अगर (args.n < 0)               info =  6;
    अगर (args.m < 0)               info =  5;
    अगर (unit < 0)                 info =  4;
    अगर (trans < 0)                info =  3;
    अगर (uplo  < 0)                info =  2;
    अगर (side  < 0)                info =  1;
  }

  अगर (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    अगर (Side == Cblas_बायाँ)         side  = 1;
    अगर (Side == Cblas_दायाँ)        side  = 0;

    अगर (Uplo == CblasUpper)        uplo  = 1;
    अगर (Uplo == CblasLower)        uplo  = 0;

    अगर (Trans == CblasNoTrans)     trans = 0;
    अगर (Trans == CblasTrans)       trans = 1;
#अगर_अघोषित COMPLEX
    अगर (Trans == CblasConjNoTrans) trans = 0;
    अगर (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    अगर (Trans == CblasConjNoTrans) trans = 2;
    अगर (Trans == CblasConjTrans)   trans = 3;
#पूर्ण_अगर

    अगर (Diag == CblasUnit)         unit  = 0;
    अगर (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    अगर (side & 1) nrowa = args.n;

    अगर (args.ldb < MAX(1,args.m)) info = 11;
    अगर (args.lda < MAX(1,nrowa))  info =  9;
    अगर (args.n < 0)               info =  6;
    अगर (args.m < 0)               info =  5;
    अगर (unit < 0)                 info =  4;
    अगर (trans < 0)                info =  3;
    अगर (uplo  < 0)                info =  2;
    अगर (side  < 0)                info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

  अगर (args.m == 0 || args.n == 0) वापस;

#अगर defined(TRMM)
#अगर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#अगर defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर
  अगर (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_बायाँ && m == lda && n == ldb) {
    अगर (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      वापस;
    } अन्यथा अगर (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      वापस;
    }
  }
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर



//end of the ifndef CBLAS ... else ... section	
#पूर्ण_अगर

  अगर ((args.m == 0) || (args.n == 0)) वापस;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  बफर = (FLOAT *)blas_memory_alloc(0);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#अगर_घोषित SMP
  mode |= (trans << BLAS_TRANSA_SHIFT);
  mode |= (side  << BLAS_RSIDE_SHIFT);

/*
  if ( args.m < 2 * GEMM_MULTITHREAD_THRESHOLD )
	args.nthreads = 1;
  else
	if ( args.n < 2 * GEMM_MULTITHREAD_THRESHOLD )
		args.nthreads = 1;
*/
  अगर ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  अन्यथा
	args.nthreads = अंक_cpu_avail(3);
		

  अगर (args.nthreads == 1) {
#पूर्ण_अगर

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP
  } अन्यथा {
    अगर (!side) {
      gemm_thread_n(mode, &args, शून्य, शून्य, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } अन्यथा {
      gemm_thread_m(mode, &args, शून्य, शून्य, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  वापस;
}

