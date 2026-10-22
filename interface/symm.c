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

#अगर_अघोषित COMPLEX
#घोषणा SMP_THRESHOLD_MIN 65536.
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QSYMM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DSYMM "
#अन्यथा
#घोषणा ERROR_NAME "SSYMM "
#पूर्ण_अगर
#अन्यथा
#घोषणा SMP_THRESHOLD_MIN 8192.
#अगर_अघोषित GEMM3M
#अगर_अघोषित HEMM
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XSYMM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZSYMM "
#अन्यथा
#घोषणा ERROR_NAME "CSYMM "
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XHEMM "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHEMM "
#अन्यथा
#घोषणा ERROR_NAME "CHEMM "
#पूर्ण_अगर
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित HEMM
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XSYMM3M "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZSYMM3M "
#अन्यथा
#घोषणा ERROR_NAME "CSYMM3M "
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XHEMM3M "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZHEMM3M "
#अन्यथा
#घोषणा ERROR_NAME "CHEMM3M "
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

#अगर_अघोषित GEMM_MULTITHREAD_THRESHOLD
#घोषणा GEMM_MULTITHREAD_THRESHOLD 4
#पूर्ण_अगर

#अगर_घोषित DYNAMIC_ARCH
बाह्य अक्षर* gotoblas_corename(व्योम);                                
#पूर्ण_अगर                                                               

#अगर_घोषित SMP
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा MODE	(BLAS_XDOUBLE | BLAS_REAL)
#या_अगर defined(DOUBLE)
#घोषणा MODE	(BLAS_DOUBLE  | BLAS_REAL)
#अन्यथा
#घोषणा MODE	(BLAS_SINGLE  | BLAS_REAL)
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा MODE	(BLAS_XDOUBLE | BLAS_COMPLEX)
#या_अगर defined(DOUBLE)
#घोषणा MODE	(BLAS_DOUBLE  | BLAS_COMPLEX)
#अन्यथा
#घोषणा MODE	(BLAS_SINGLE  | BLAS_COMPLEX)
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर

अटल पूर्णांक (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#अगर_अघोषित GEMM3M
#अगर_अघोषित HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#पूर्ण_अगर
#अन्यथा
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#पूर्ण_अगर
#पूर्ण_अगर
#अन्यथा
#अगर_अघोषित HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#पूर्ण_अगर
#अन्यथा
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#अगर defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
};

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *SIDE, अक्षर *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  अक्षर side_arg  = *SIDE;
  अक्षर uplo_arg  = *UPLO;

  blas_arg_t args;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर defined(SMP) && !defined(NO_AFFINITY)
  पूर्णांक nodes;
#पूर्ण_अगर
# अगर defined(SMP)
  द्विगुन MN;
#पूर्ण_अगर
  blasint info;
  पूर्णांक side;
  पूर्णांक uplo;

  PRINT_DEBUG_NAME;

  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  अगर (side_arg  == 'L') side  = 0;
  अगर (side_arg  == 'R') side  = 1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;

  args.m = *M;
  args.n = *N;

  args.c = (व्योम *)c;
  args.ldc = *ldC;

  info = 0;

  अगर (args.ldc < MAX(1, args.m)) info = 12;

  अगर (!side) {
    args.a = (व्योम *)a;
    args.b = (व्योम *)b;

    args.lda = *ldA;
    args.ldb = *ldB;

    अगर (args.ldb < MAX(1, args.m)) info =  9;
    अगर (args.lda < MAX(1, args.m)) info =  7;

  } अन्यथा {
    args.a = (व्योम *)b;
    args.b = (व्योम *)a;

    args.lda = *ldB;
    args.ldb = *ldA;

  अगर (args.lda < MAX(1, args.m)) info =  9;
  अगर (args.ldb < MAX(1, args.n)) info =  7;
  }

  अगर (args.n   < 0)              info =  4;
  अगर (args.m   < 0)              info =  3;
  अगर (uplo     < 0)              info =  2;
  अगर (side     < 0)              info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_SIDE Side, क्रमागत CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   व्योम *valpha,
	   व्योम *va, blasint lda,
	   व्योम *vb, blasint ldb,
	   व्योम *vbeta,
	   व्योम *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#पूर्ण_अगर

  blas_arg_t args;
  पूर्णांक side, uplo;
  blasint info;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#अगर defined(SMP) && !defined(NO_AFFINITY)
  पूर्णांक nodes;
#पूर्ण_अगर
#अगर defined(SMP)
  द्विगुन MN;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

#अगर_अघोषित COMPLEX
  args.alpha = (व्योम *)&alpha;
  args.beta  = (व्योम *)&beta;
#अन्यथा
  args.alpha = (व्योम *)alpha;
  args.beta  = (व्योम *)beta;
#पूर्ण_अगर

  args.c = (व्योम *)c;
  args.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  अगर (order == CblasColMajor) {
    अगर (Side == Cblas_बायाँ)  side = 0;
    अगर (Side == Cblas_दायाँ) side = 1;

    अगर (Uplo == CblasUpper) uplo  = 0;
    अगर (Uplo == CblasLower) uplo  = 1;

    info = -1;

    args.m = m;
    args.n = n;

    अगर (args.ldc < MAX(1, args.m)) info = 12;

    अगर (!side) {
      args.a = (व्योम *)a;
      args.b = (व्योम *)b;

      args.lda = lda;
      args.ldb = ldb;

      अगर (args.ldb < MAX(1, args.m)) info =  9;
      अगर (args.lda < MAX(1, args.m)) info =  7;

    } अन्यथा {
      args.a = (व्योम *)b;
      args.b = (व्योम *)a;

      args.lda = ldb;
      args.ldb = lda;

      अगर (args.lda < MAX(1, args.m)) info =  9;
      अगर (args.ldb < MAX(1, args.n)) info =  7;
    }

    अगर (args.n   < 0)              info =  4;
    अगर (args.m   < 0)              info =  3;
    अगर (uplo     < 0)              info =  2;
    अगर (side     < 0)              info =  1;
  }

  अगर (order == CblasRowMajor) {
    अगर (Side == Cblas_बायाँ)  side = 1;
    अगर (Side == Cblas_दायाँ) side = 0;

    अगर (Uplo == CblasUpper) uplo  = 1;
    अगर (Uplo == CblasLower) uplo  = 0;

    info = -1;

    args.m = n;
    args.n = m;

    अगर (args.ldc < MAX(1, args.m)) info = 12;

    अगर (!side) {
      args.a = (व्योम *)a;
      args.b = (व्योम *)b;

      args.lda = lda;
      args.ldb = ldb;

      अगर (args.ldb < MAX(1, args.m)) info =  9;
      अगर (args.lda < MAX(1, args.m)) info =  7;

    } अन्यथा {
      args.a = (व्योम *)b;
      args.b = (व्योम *)a;

      args.lda = ldb;
      args.ldb = lda;

      अगर (args.lda < MAX(1, args.m)) info =  9;
      अगर (args.ldb < MAX(1, args.n)) info =  7;
    }

    अगर (args.n   < 0)              info =  4;
    अगर (args.m   < 0)              info =  3;
    अगर (uplo     < 0)              info =  2;
    अगर (side     < 0)              info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

   अगर (args.m == 0 || args.n == 0) वापस;
#अगर !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#अगर defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#अगर defined(DYNAMIC_ARCH)
अगर (म_भेद(gotoblas_corename(), "armv9sme") == 0
#अगर defined(__clang__)
 || म_भेद(gotoblas_corename(), "vortexm4") == 0
#पूर्ण_अगर
)
#पूर्ण_अगर
   अगर (order == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     अगर (Side == Cblas_बायाँ && Uplo == CblasUpper) {
       SSYMM_DIRECT_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); वापस;
     }
     अन्यथा अगर (Side == Cblas_बायाँ && Uplo == CblasLower) {
       SSYMM_DIRECT_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); वापस;
     }
   }
#पूर्ण_अगर
#पूर्ण_अगर

#पूर्ण_अगर

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
  args.common = शून्य;
  MN = 2.* (द्विगुन) args.m * (द्विगुन)args.m * (द्विगुन) args.n;
  अगर (MN <= (SMP_THRESHOLD_MIN * (द्विगुन) GEMM_MULTITHREAD_THRESHOLD) ) {
	  args.nthreads = 1;
  } अन्यथा {
	 args.nthreads = अंक_cpu_avail(3);
  }
  अगर (args.nthreads == 1) {
#पूर्ण_अगर

    (symm[(side << 1) | uplo ])(&args, शून्य, शून्य, sa, sb, 0);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_अघोषित NO_AFFINITY
    nodes = get_अंक_nodes();

    अगर (nodes > 1) {

      args.nthreads /= nodes;

      gemm_thread_mn(MODE, &args, शून्य, शून्य,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } अन्यथा {
#पूर्ण_अगर

#अगर_अघोषित USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&args, शून्य, शून्य, sa, sb, 0);

#अन्यथा

      GEMM_THREAD(MODE, &args, शून्य, शून्य, symm[(side << 1) | uplo ], sa, sb, args.nthreads);

#पूर्ण_अगर

#अगर_अघोषित NO_AFFINITY
    }
#पूर्ण_अगर

  }
#पूर्ण_अगर

 blas_memory_free(बफर);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE,
		       (!side)? args.m * (args.m / 2 + args.n) : args.n * (args.m + args.n / 2),
		       (!side)? 2 * args.m * args.m * args.n : 2 * args.m * args.n * args.n);

  IDEBUG_END;

  वापस;
}
