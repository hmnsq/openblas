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

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#अन्यथा
#define ERROR_NAME "SSYRK "
#endif
#अन्यथा
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#अन्यथा
#define ERROR_NAME "CSYRK "
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#अन्यथा
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

स्थिर पूर्णाङ्क (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#अन्यथा
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  द्विगुण NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क trans;
  पूर्णाङ्क nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (शून्य *)a;
  args.c = (शून्य *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
#ifdef HEMM
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'C') trans = 1;
#अन्यथा
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  यदि (trans & 1) nrowa = args.k;

  info = 0;

  यदि (args.ldc < MAX(1,args.n)) info = 10;
  यदि (args.lda < MAX(1,nrowa))  info =  7;
  यदि (args.k < 0)               info =  4;
  यदि (args.n < 0)               info =  3;
  यदि (trans < 0)                info =  2;
  यदि (uplo  < 0)                info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, गणना CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#यदि !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#अन्यथा
	   शून्य *valpha,
#endif
#यदि !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#अन्यथा
	   शून्य *va, blasint lda,
#endif
#यदि !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#अन्यथा
	   शून्य *vbeta,
#endif
#यदि !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   शून्य *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#यदि !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  पूर्णाङ्क uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
द्विगुण NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  पूर्णाङ्क mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  पूर्णाङ्क mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  पूर्णाङ्क mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (शून्य *)a;
  args.c = (शून्य *)c;

  args.lda = lda;
  args.ldc = ldc;

#यदि !defined(COMPLEX) || defined(HEMM)
  args.alpha = (शून्य *)&alpha;
  args.beta  = (शून्य *)&beta;
#अन्यथा
  args.alpha = (शून्य *)alpha;
  args.beta  = (शून्य *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  यदि (order == CblasColMajor) {
    यदि (Uplo == CblasUpper) uplo  = 0;
    यदि (Uplo == CblasLower) uplo  = 1;

    यदि (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    यदि (Trans == CblasTrans)       trans = 1;
    यदि (Trans == CblasConjNoTrans) trans = 0;
    यदि (Trans == CblasConjTrans)   trans = 1;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 1;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    यदि (trans & 1) nrowa = args.k;

    यदि (args.ldc < MAX(1,args.n)) info = 10;
    यदि (args.lda < MAX(1,nrowa))  info =  7;
    यदि (args.k < 0)               info =  4;
    यदि (args.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (order == CblasRowMajor) {
    यदि (Uplo == CblasUpper) uplo  = 1;
    यदि (Uplo == CblasLower) uplo  = 0;

    यदि (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    यदि (Trans == CblasTrans)       trans = 0;
    यदि (Trans == CblasConjNoTrans) trans = 1;
    यदि (Trans == CblasConjTrans)   trans = 0;
#नोचेद्यदि !defined(HEMM)
    यदि (Trans == CblasTrans)       trans = 0;
#अन्यथा
    यदि (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    यदि (trans & 1) nrowa = args.k;

    यदि (args.ldc < MAX(1,args.n)) info = 10;
    यदि (args.lda < MAX(1,nrowa))  info =  7;
    यदि (args.k < 0)               info =  4;
    यदि (args.n < 0)               info =  3;
    यदि (trans < 0)                info =  2;
    यदि (uplo  < 0)                info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (args.n == 0) प्रत्यागम;

#ifdef DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);
#endif
	
#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (म_भेद(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  यदि (order == CblasRowMajor && n == ldc) {
    यदि (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      प्रत्यागम;
    } अन्यथा यदि (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      प्रत्यागम;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  यदि (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  यदि (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } अन्यथा {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = रिक्त;

  NNK = (द्विगुण)(args.n+1)*(द्विगुण)args.n*(द्विगुण)args.k;
  यदि (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } अन्यथा {
  args.nthreads = अंक_cpu_avail(3);
  }

  यदि (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, रिक्त, रिक्त, sa, sb, 0);

#ifdef SMP
  } अन्यथा {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

    syrk_thread(mode, &args, रिक्त, रिक्त, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  प्रत्यागम;
}
