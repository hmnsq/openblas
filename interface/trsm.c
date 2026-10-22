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
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifndef TRMM
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRSM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DTRSM "
#अन्यथा
#define ERROR_NAME "STRSM "
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XTRSM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZTRSM "
#अन्यथा
#define ERROR_NAME "CTRSM "
#endif
#endif
#अन्यथा
#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QTRMM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DTRMM "
#अन्यथा
#define ERROR_NAME "STRMM "
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XTRMM "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZTRMM "
#अन्यथा
#define ERROR_NAME "CTRMM "
#endif
#endif
#endif

#ifndef COMPLEX
#define SMP_FACTOR 256
#अन्यथा
#define SMP_FACTOR 128
#endif

#ifdef DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(शून्य);                                
#endif                                                               

स्थिर पूर्णांक (*trsm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef TRMM
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
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *SIDE, अक्षर *UPLO, अक्षर *TRANS, अक्षर *DIAG,
	   blasint *M, blasint *N, FLOAT *alpha,
	   FLOAT *a, blasint *ldA, FLOAT *b, blasint *ldB){

  अक्षर side_arg  = *SIDE;
  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;
  अक्षर diag_arg  = *DIAG;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  blasint info;
  पूर्णांक side;
  पूर्णांक uplo;
  पूर्णांक unit;
  पूर्णांक trans;
  पूर्णांक nrowa;

  PRINT_DEBUG_NAME;

  args.m = *M;
  args.n = *N;

  args.a = (शून्य *)a;
  args.b = (शून्य *)b;

  args.lda = *ldA;
  args.ldb = *ldB;

  args.beta = (शून्य *)alpha;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  side  = -1;
  trans = -1;
  unit  = -1;
  uplo  = -1;

  जर (side_arg  == 'L') side  = 0;
  जर (side_arg  == 'R') side  = 1;

  जर (trans_arg == 'N') trans = 0;
  जर (trans_arg == 'T') trans = 1;
  जर (trans_arg == 'R') trans = 2;
  जर (trans_arg == 'C') trans = 3;

  जर (diag_arg  == 'U') unit  = 0;
  जर (diag_arg  == 'N') unit  = 1;

  जर (uplo_arg  == 'U') uplo  = 0;
  जर (uplo_arg  == 'L') uplo  = 1;

  nrowa = args.m;
  जर (side & 1) nrowa = args.n;

  info = 0;

  जर (args.ldb < MAX(1,args.m)) info = 11;
  जर (args.lda < MAX(1,nrowa))  info =  9;
  जर (args.n < 0)               info =  6;
  जर (args.m < 0)               info =  5;
  जर (unit < 0)                 info =  4;
  जर (trans < 0)                info =  3;
  जर (uplo  < 0)                info =  2;
  जर (side  < 0)                info =  1;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME)-1);
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_SIDE Side,  गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE Trans, गणना CBLAS_DIAG Diag,
	   blasint m, blasint n,
#ifndef COMPLEX
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
#endif

  blas_arg_t args;
  पूर्णांक side, uplo, trans, unit;
  blasint info, nrowa;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णांक mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  पूर्णांक mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णांक mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.a = (शून्य *)a;
  args.b = (शून्य *)b;

  args.lda = lda;
  args.ldb = ldb;

#ifndef COMPLEX
  args.beta = (शून्य *)&alpha;
#अन्यथा
  args.beta = (शून्य *)alpha;
#endif

  side   = -1;
  uplo   = -1;
  trans  = -1;
  unit   = -1;
  info   =  0;

  जर (order == CblasColMajor) {
    args.m = m;
    args.n = n;

    जर (Side == Cblas_डावा)         side  = 0;
    जर (Side == Cblas_उजवा)        side  = 1;

    जर (Uplo == CblasUpper)        uplo  = 0;
    जर (Uplo == CblasLower)        uplo  = 1;

    जर (Trans == CblasNoTrans)     trans = 0;
    जर (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    जर (Trans == CblasConjNoTrans) trans = 0;
    जर (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    जर (Trans == CblasConjNoTrans) trans = 2;
    जर (Trans == CblasConjTrans)   trans = 3;
#endif

    जर (Diag == CblasUnit)          unit  = 0;
    जर (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    nrowa = args.m;
    जर (side & 1) nrowa = args.n;

    जर (args.ldb < MAX(1,args.m)) info = 11;
    जर (args.lda < MAX(1,nrowa))  info =  9;
    जर (args.n < 0)               info =  6;
    जर (args.m < 0)               info =  5;
    जर (unit < 0)                 info =  4;
    जर (trans < 0)                info =  3;
    जर (uplo  < 0)                info =  2;
    जर (side  < 0)                info =  1;
  }

  जर (order == CblasRowMajor) {
    args.m = n;
    args.n = m;

    जर (Side == Cblas_डावा)         side  = 1;
    जर (Side == Cblas_उजवा)        side  = 0;

    जर (Uplo == CblasUpper)        uplo  = 1;
    जर (Uplo == CblasLower)        uplo  = 0;

    जर (Trans == CblasNoTrans)     trans = 0;
    जर (Trans == CblasTrans)       trans = 1;
#ifndef COMPLEX
    जर (Trans == CblasConjNoTrans) trans = 0;
    जर (Trans == CblasConjTrans)   trans = 1;
#अन्यथा
    जर (Trans == CblasConjNoTrans) trans = 2;
    जर (Trans == CblasConjTrans)   trans = 3;
#endif

    जर (Diag == CblasUnit)         unit  = 0;
    जर (Diag == CblasNonUnit)      unit  = 1;

    info = -1;

    nrowa = args.m;
    जर (side & 1) nrowa = args.n;

    जर (args.ldb < MAX(1,args.m)) info = 11;
    जर (args.lda < MAX(1,nrowa))  info =  9;
    जर (args.n < 0)               info =  6;
    जर (args.m < 0)               info =  5;
    जर (unit < 0)                 info =  4;
    जर (trans < 0)                info =  3;
    जर (uplo  < 0)                info =  2;
    जर (side  < 0)                info =  1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

  जर (args.m == 0 || args.n == 0) परत;

#जर defined(TRMM)
#जर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#जर defined(ARCH_ARM64) && (defined(USE_STRMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#जर defined(DYNAMIC_ARCH)
जर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#जर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  जर (order == CblasRowMajor && Diag == CblasNonUnit && Side == Cblas_डावा && m == lda && n == ldb) {
    जर (Trans == CblasNoTrans || Trans == CblasConjNoTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LNUN : STRMM_DIRECT_LNLN)(m, n, alpha, a, lda, b, ldb);
      परत;
    } अन्यथा जर (Trans == CblasTrans || Trans == CblasConjTrans) {
      (Uplo == CblasUpper ? STRMM_DIRECT_LTUN : STRMM_DIRECT_LTLN)(m, n, alpha, a, lda, b, ldb);
      परत;
    }
  }
#endif
#endif
#endif



//end of the ifndef CBLAS ... else ... section	
#endif

  जर ((args.m == 0) || (args.n == 0)) परत;
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
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
  जर ( args.m * args.n < SMP_FACTOR * GEMM_MULTITHREAD_THRESHOLD)
	args.nthreads = 1;
  अन्यथा
	args.nthreads = अंक_cpu_avail(3);
		

  जर (args.nthreads == 1) {
#endif

    (trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP
  } अन्यथा {
    जर (!side) {
      gemm_thread_n(mode, &args, रिक्त, रिक्त, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    } अन्यथा {
      gemm_thread_m(mode, &args, रिक्त, रिक्त, trsm[(side<<4) | (trans<<2) | (uplo<<1) | unit], sa, sb, args.nthreads);
    }
  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side) ? args.m * (args.m + args.n) : args.n * (args.m + args.n),
		       (!side) ? args.m * args.m * args.n : args.m * args.n * args.n);

  IDEBUG_END;

  परत;
}

