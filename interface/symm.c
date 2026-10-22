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

#यदि_अपरिभाषितम् COMPLEX
#परिभाषय SMP_देहली_MIN 65536.
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QSYMM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DSYMM "
#अन्यथा
#परिभाषय ERROR_NAME "SSYMM "
#यदि_अन्तः
#अन्यथा
#परिभाषय SMP_देहली_MIN 8192.
#यदि_अपरिभाषितम् GEMM3M
#यदि_अपरिभाषितम् HEMM
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XSYMM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZSYMM "
#अन्यथा
#परिभाषय ERROR_NAME "CSYMM "
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XHEMM "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZHEMM "
#अन्यथा
#परिभाषय ERROR_NAME "CHEMM "
#यदि_अन्तः
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् HEMM
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XSYMM3M "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZSYMM3M "
#अन्यथा
#परिभाषय ERROR_NAME "CSYMM3M "
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XHEMM3M "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZHEMM3M "
#अन्यथा
#परिभाषय ERROR_NAME "CHEMM3M "
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

#यदि_अपरिभाषितम् GEMM_MULTITHREAD_देहली
#परिभाषय GEMM_MULTITHREAD_देहली 4
#यदि_अन्तः

#यदि_परिभाषितम् DYNAMIC_ARCH
बाह्य वर्ण* gotoblas_corename(शून्य);                                
#यदि_अन्तः                                                               

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय बहुलकम्_2	(BLAS_XDOUBLE | BLAS_REAL)
#नोचेद्यदि defined(DOUBLE)
#परिभाषय बहुलकम्_2	(BLAS_DOUBLE  | BLAS_REAL)
#अन्यथा
#परिभाषय बहुलकम्_2	(BLAS_SINGLE  | BLAS_REAL)
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय बहुलकम्_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#नोचेद्यदि defined(DOUBLE)
#परिभाषय बहुलकम्_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#अन्यथा
#परिभाषय बहुलकम्_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः

स्थिर पूर्णाङ्क (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#यदि_अपरिभाषितम् GEMM3M
#यदि_अपरिभाषितम् HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#यदि_अन्तः
#अन्यथा
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#यदि_अन्तः
#यदि_अन्तः
#अन्यथा
#यदि_अपरिभाषितम् HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#यदि_अन्तः
#अन्यथा
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#यदि defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
};

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *SIDE, वर्ण *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  वर्ण side_arg  = *SIDE;
  वर्ण uplo_arg  = *UPLO;

  blas_arg_t तर्काः;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि defined(SMP) && !defined(NO_AFFINITY)
  पूर्णाङ्क nodes;
#यदि_अन्तः
# यदि defined(SMP)
  द्विगुण MN;
#यदि_अन्तः
  blasint info;
  पूर्णाङ्क side;
  पूर्णाङ्क uplo;

  PRINT_DEBUG_NAME;

  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  यदि (side_arg  == 'L') side  = 0;
  यदि (side_arg  == 'R') side  = 1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  तर्काः.m = *M;
  तर्काः.n = *N;

  तर्काः.c = (शून्य *)c;
  तर्काः.ldc = *ldC;

  info = 0;

  यदि (तर्काः.ldc < MAX(1, तर्काः.m)) info = 12;

  यदि (!side) {
    तर्काः.a = (शून्य *)a;
    तर्काः.b = (शून्य *)b;

    तर्काः.lda = *ldA;
    तर्काः.ldb = *ldB;

    यदि (तर्काः.ldb < MAX(1, तर्काः.m)) info =  9;
    यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  7;

  } अन्यथा {
    तर्काः.a = (शून्य *)b;
    तर्काः.b = (शून्य *)a;

    तर्काः.lda = *ldB;
    तर्काः.ldb = *ldA;

  यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  9;
  यदि (तर्काः.ldb < MAX(1, तर्काः.n)) info =  7;
  }

  यदि (तर्काः.n   < 0)              info =  4;
  यदि (तर्काः.m   < 0)              info =  3;
  यदि (uplo     < 0)              info =  2;
  यदि (side     < 0)              info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_SIDE Side, गणना CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#अन्यथा
	   शून्य *valpha,
	   शून्य *va, blasint lda,
	   शून्य *vb, blasint ldb,
	   शून्य *vbeta,
	   शून्य *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#यदि_अन्तः

  blas_arg_t तर्काः;
  पूर्णाङ्क side, uplo;
  blasint info;

  FLOAT *बफर;
  FLOAT *sa, *sb;

#यदि defined(SMP) && !defined(NO_AFFINITY)
  पूर्णाङ्क nodes;
#यदि_अन्तः
#यदि defined(SMP)
  द्विगुण MN;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

#यदि_अपरिभाषितम् COMPLEX
  तर्काः.alpha = (शून्य *)&alpha;
  तर्काः.beta  = (शून्य *)&beta;
#अन्यथा
  तर्काः.alpha = (शून्य *)alpha;
  तर्काः.beta  = (शून्य *)beta;
#यदि_अन्तः

  तर्काः.c = (शून्य *)c;
  तर्काः.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (Side == Cblas_वामः)  side = 0;
    यदि (Side == Cblas_दक्षिणः) side = 1;

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 1;

    info = -1;

    तर्काः.m = m;
    तर्काः.n = n;

    यदि (तर्काः.ldc < MAX(1, तर्काः.m)) info = 12;

    यदि (!side) {
      तर्काः.a = (शून्य *)a;
      तर्काः.b = (शून्य *)b;

      तर्काः.lda = lda;
      तर्काः.ldb = ldb;

      यदि (तर्काः.ldb < MAX(1, तर्काः.m)) info =  9;
      यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  7;

    } अन्यथा {
      तर्काः.a = (शून्य *)b;
      तर्काः.b = (शून्य *)a;

      तर्काः.lda = ldb;
      तर्काः.ldb = lda;

      यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  9;
      यदि (तर्काः.ldb < MAX(1, तर्काः.n)) info =  7;
    }

    यदि (तर्काः.n   < 0)              info =  4;
    यदि (तर्काः.m   < 0)              info =  3;
    यदि (uplo     < 0)              info =  2;
    यदि (side     < 0)              info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (Side == Cblas_वामः)  side = 1;
    यदि (Side == Cblas_दक्षिणः) side = 0;

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 0;

    info = -1;

    तर्काः.m = n;
    तर्काः.n = m;

    यदि (तर्काः.ldc < MAX(1, तर्काः.m)) info = 12;

    यदि (!side) {
      तर्काः.a = (शून्य *)a;
      तर्काः.b = (शून्य *)b;

      तर्काः.lda = lda;
      तर्काः.ldb = ldb;

      यदि (तर्काः.ldb < MAX(1, तर्काः.m)) info =  9;
      यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  7;

    } अन्यथा {
      तर्काः.a = (शून्य *)b;
      तर्काः.b = (शून्य *)a;

      तर्काः.lda = ldb;
      तर्काः.ldb = lda;

      यदि (तर्काः.lda < MAX(1, तर्काः.m)) info =  9;
      यदि (तर्काः.ldb < MAX(1, तर्काः.n)) info =  7;
    }

    यदि (तर्काः.n   < 0)              info =  4;
    यदि (तर्काः.m   < 0)              info =  3;
    यदि (uplo     < 0)              info =  2;
    यदि (side     < 0)              info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

   यदि (तर्काः.m == 0 || तर्काः.n == 0) प्रत्यागम;
#यदि !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#यदि defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_प्रत्यक्षम्)||defined(DYNAMIC_ARCH))
#यदि defined(DYNAMIC_ARCH)
यदि (वाक्यं_तुलय(gotoblas_corename(), "armv9sme") == 0
#यदि defined(__clang__)
 || वाक्यं_तुलय(gotoblas_corename(), "vortexm4") == 0
#यदि_अन्तः
)
#यदि_अन्तः
   यदि (क्रमः == CblasRowMajor && m == lda && n == ldb && n == ldc)
   {
     यदि (Side == Cblas_वामः && Uplo == Cblas_ऊर्ध्वस्थम्) {
       SSYMM_प्रत्यक्षम्_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); प्रत्यागम;
     }
     अन्यथा यदि (Side == Cblas_वामः && Uplo == Cblas_अधस्तनम्) {
       SSYMM_प्रत्यक्षम्_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); प्रत्यागम;
     }
   }
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः

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
  तर्काः.common = रिक्त;
  MN = 2.* (द्विगुण) तर्काः.m * (द्विगुण)तर्काः.m * (द्विगुण) तर्काः.n;
  यदि (MN <= (SMP_देहली_MIN * (द्विगुण) GEMM_MULTITHREAD_देहली) ) {
	  तर्काः.nthreads = 1;
  } अन्यथा {
	 तर्काः.nthreads = अङ्कः_cpu_avail(3);
  }
  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    (symm[(side << 1) | uplo ])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_अपरिभाषितम् NO_AFFINITY
    nodes = get_अङ्कः_nodes();

    यदि (nodes > 1) {

      तर्काः.nthreads /= nodes;

      gemm_thread_mn(बहुलकम्_2, &तर्काः, रिक्त, रिक्त,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } अन्यथा {
#यदि_अन्तः

#यदि_अपरिभाषितम् USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#अन्यथा

      GEMM_THREAD(बहुलकम्_2, &तर्काः, रिक्त, रिक्त, symm[(side << 1) | uplo ], sa, sb, तर्काः.nthreads);

#यदि_अन्तः

#यदि_अपरिभाषितम् NO_AFFINITY
    }
#यदि_अन्तः

  }
#यदि_अन्तः

 blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE,
		       (!side)? तर्काः.m * (तर्काः.m / 2 + तर्काः.n) : तर्काः.n * (तर्काः.m + तर्काः.n / 2),
		       (!side)? 2 * तर्काः.m * तर्काः.m * तर्काः.n : 2 * तर्काः.m * तर्काः.n * तर्काः.n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
