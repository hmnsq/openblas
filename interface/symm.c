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

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచించు SMP_పరిమితి_MIN 65536.
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QSYMM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DSYMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "SSYMM "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచించు SMP_పరిమితి_MIN 8192.
#నిర్వచితం_కాకపోతే GEMM3M
#నిర్వచితం_కాకపోతే HEMM
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XSYMM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZSYMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CSYMM "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XHEMM "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZHEMM "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CHEMM "
#అయితే_ముగింపు
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితం_కాకపోతే HEMM
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XSYMM3M "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZSYMM3M "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CSYMM3M "
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XHEMM3M "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZHEMM3M "
#లేకపోతే
#నిర్వచించు ERROR_NAME "CHEMM3M "
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే GEMM_MULTITHREAD_పరిమితి
#నిర్వచించు GEMM_MULTITHREAD_పరిమితి 4
#అయితే_ముగింపు

#నిర్వచితమైతే DYNAMIC_ARCH
బాహ్య అక్షరం* gotoblas_corename(శూన్య);                                
#అయితే_ముగింపు                                                               

#నిర్వచితమైతే SMP
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు బహుళకం_2	(BLAS_XDOUBLE | BLAS_REAL)
#లేదావేళ defined(DOUBLE)
#నిర్వచించు బహుళకం_2	(BLAS_DOUBLE  | BLAS_REAL)
#లేకపోతే
#నిర్వచించు బహుళకం_2	(BLAS_SINGLE  | BLAS_REAL)
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు బహుళకం_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#లేదావేళ defined(DOUBLE)
#నిర్వచించు బహుళకం_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#లేకపోతే
#నిర్వచించు బహుళకం_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు

స్థిర పూర్ణాంకం (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#నిర్వచితం_కాకపోతే GEMM3M
#నిర్వచితం_కాకపోతే HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#అయితే_ముగింపు
#లేకపోతే
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#అయితే_ముగింపు
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితం_కాకపోతే HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#అయితే_ముగింపు
#లేకపోతే
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#ఒకవేళ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
};

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(అక్షరం *SIDE, అక్షరం *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  అక్షరం side_arg  = *SIDE;
  అక్షరం uplo_arg  = *UPLO;

  blas_arg_t వాదనలు;

  FLOAT *బఫరు;
  FLOAT *sa, *sb;

#ఒకవేళ defined(SMP) && !defined(NO_AFFINITY)
  పూర్ణాంకం nodes;
#అయితే_ముగింపు
# ఒకవేళ defined(SMP)
  ద్విగుణ MN;
#అయితే_ముగింపు
  blasint info;
  పూర్ణాంకం side;
  పూర్ణాంకం uplo;

  PRINT_DEBUG_NAME;

  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  ఒకవేళ (side_arg  == 'L') side  = 0;
  ఒకవేళ (side_arg  == 'R') side  = 1;

  ఒకవేళ (uplo_arg  == 'U') uplo  = 0;
  ఒకవేళ (uplo_arg  == 'L') uplo  = 1;

  వాదనలు.m = *M;
  వాదనలు.n = *N;

  వాదనలు.c = (శూన్య *)c;
  వాదనలు.ldc = *ldC;

  info = 0;

  ఒకవేళ (వాదనలు.ldc < MAX(1, వాదనలు.m)) info = 12;

  ఒకవేళ (!side) {
    వాదనలు.a = (శూన్య *)a;
    వాదనలు.b = (శూన్య *)b;

    వాదనలు.lda = *ldA;
    వాదనలు.ldb = *ldB;

    ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.m)) info =  9;
    ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  7;

  } లేకపోతే {
    వాదనలు.a = (శూన్య *)b;
    వాదనలు.b = (శూన్య *)a;

    వాదనలు.lda = *ldB;
    వాదనలు.ldb = *ldA;

  ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  9;
  ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.n)) info =  7;
  }

  ఒకవేళ (వాదనలు.n   < 0)              info =  4;
  ఒకవేళ (వాదనలు.m   < 0)              info =  3;
  ఒకవేళ (uplo     < 0)              info =  2;
  ఒకవేళ (side     < 0)              info =  1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_SIDE Side, గణన CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#నిర్వచితం_కాకపోతే COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#లేకపోతే
	   శూన్య *valpha,
	   శూన్య *va, blasint lda,
	   శూన్య *vb, blasint ldb,
	   శూన్య *vbeta,
	   శూన్య *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#అయితే_ముగింపు

  blas_arg_t వాదనలు;
  పూర్ణాంకం side, uplo;
  blasint info;

  FLOAT *బఫరు;
  FLOAT *sa, *sb;

#ఒకవేళ defined(SMP) && !defined(NO_AFFINITY)
  పూర్ణాంకం nodes;
#అయితే_ముగింపు
#ఒకవేళ defined(SMP)
  ద్విగుణ MN;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

#నిర్వచితం_కాకపోతే COMPLEX
  వాదనలు.alpha = (శూన్య *)&alpha;
  వాదనలు.beta  = (శూన్య *)&beta;
#లేకపోతే
  వాదనలు.alpha = (శూన్య *)alpha;
  వాదనలు.beta  = (శూన్య *)beta;
#అయితే_ముగింపు

  వాదనలు.c = (శూన్య *)c;
  వాదనలు.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
    ఒకవేళ (Side == Cblas_ఎడమ)  side = 0;
    ఒకవేళ (Side == Cblas_కుడి) side = 1;

    ఒకవేళ (Uplo == Cblas_పై) uplo  = 0;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 1;

    info = -1;

    వాదనలు.m = m;
    వాదనలు.n = n;

    ఒకవేళ (వాదనలు.ldc < MAX(1, వాదనలు.m)) info = 12;

    ఒకవేళ (!side) {
      వాదనలు.a = (శూన్య *)a;
      వాదనలు.b = (శూన్య *)b;

      వాదనలు.lda = lda;
      వాదనలు.ldb = ldb;

      ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.m)) info =  9;
      ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  7;

    } లేకపోతే {
      వాదనలు.a = (శూన్య *)b;
      వాదనలు.b = (శూన్య *)a;

      వాదనలు.lda = ldb;
      వాదనలు.ldb = lda;

      ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  9;
      ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.n)) info =  7;
    }

    ఒకవేళ (వాదనలు.n   < 0)              info =  4;
    ఒకవేళ (వాదనలు.m   < 0)              info =  3;
    ఒకవేళ (uplo     < 0)              info =  2;
    ఒకవేళ (side     < 0)              info =  1;
  }

  ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
    ఒకవేళ (Side == Cblas_ఎడమ)  side = 1;
    ఒకవేళ (Side == Cblas_కుడి) side = 0;

    ఒకవేళ (Uplo == Cblas_పై) uplo  = 1;
    ఒకవేళ (Uplo == Cblas_క్రింది) uplo  = 0;

    info = -1;

    వాదనలు.m = n;
    వాదనలు.n = m;

    ఒకవేళ (వాదనలు.ldc < MAX(1, వాదనలు.m)) info = 12;

    ఒకవేళ (!side) {
      వాదనలు.a = (శూన్య *)a;
      వాదనలు.b = (శూన్య *)b;

      వాదనలు.lda = lda;
      వాదనలు.ldb = ldb;

      ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.m)) info =  9;
      ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  7;

    } లేకపోతే {
      వాదనలు.a = (శూన్య *)b;
      వాదనలు.b = (శూన్య *)a;

      వాదనలు.lda = ldb;
      వాదనలు.ldb = lda;

      ఒకవేళ (వాదనలు.lda < MAX(1, వాదనలు.m)) info =  9;
      ఒకవేళ (వాదనలు.ldb < MAX(1, వాదనలు.n)) info =  7;
    }

    ఒకవేళ (వాదనలు.n   < 0)              info =  4;
    ఒకవేళ (వాదనలు.m   < 0)              info =  3;
    ఒకవేళ (uplo     < 0)              info =  2;
    ఒకవేళ (side     < 0)              info =  1;
  }

  ఒకవేళ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
    తిరిగి;
  }

   ఒకవేళ (వాదనలు.m == 0 || వాదనలు.n == 0) తిరిగి;
#ఒకవేళ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ఒకవేళ defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_ప్రత్యక్ష)||defined(DYNAMIC_ARCH))
#ఒకవేళ defined(DYNAMIC_ARCH)
ఒకవేళ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ఒకవేళ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#అయితే_ముగింపు
)
#అయితే_ముగింపు
   ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major && m == lda && n == ldb && n == ldc)
   {
     ఒకవేళ (Side == Cblas_ఎడమ && Uplo == Cblas_పై) {
       SSYMM_ప్రత్యక్ష_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); తిరిగి;
     }
     లేకపోతే ఒకవేళ (Side == Cblas_ఎడమ && Uplo == Cblas_క్రింది) {
       SSYMM_ప్రత్యక్ష_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); తిరిగి;
     }
   }
#అయితే_ముగింపు
#అయితే_ముగింపు

#అయితే_ముగింపు

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
  వాదనలు.common = ఖాళీ;
  MN = 2.* (ద్విగుణ) వాదనలు.m * (ద్విగుణ)వాదనలు.m * (ద్విగుణ) వాదనలు.n;
  ఒకవేళ (MN <= (SMP_పరిమితి_MIN * (ద్విగుణ) GEMM_MULTITHREAD_పరిమితి) ) {
	  వాదనలు.nthreads = 1;
  } లేకపోతే {
	 వాదనలు.nthreads = అంకె_cpu_avail(3);
  }
  ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    (symm[(side << 1) | uplo ])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP
  } లేకపోతే {

#నిర్వచితం_కాకపోతే NO_AFFINITY
    nodes = get_అంకె_nodes();

    ఒకవేళ (nodes > 1) {

      వాదనలు.nthreads /= nodes;

      gemm_thread_mn(బహుళకం_2, &వాదనలు, ఖాళీ, ఖాళీ,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } లేకపోతే {
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#లేకపోతే

      GEMM_THREAD(బహుళకం_2, &వాదనలు, ఖాళీ, ఖాళీ, symm[(side << 1) | uplo ], sa, sb, వాదనలు.nthreads);

#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే NO_AFFINITY
    }
#అయితే_ముగింపు

  }
#అయితే_ముగింపు

 blas_జ్ఞాపకం_free(బఫరు);

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE,
		       (!side)? వాదనలు.m * (వాదనలు.m / 2 + వాదనలు.n) : వాదనలు.n * (వాదనలు.m + వాదనలు.n / 2),
		       (!side)? 2 * వాదనలు.m * వాదనలు.m * వాదనలు.n : 2 * వాదనలు.m * వాదనలు.n * వాదనలు.n);

  IDEBUG_చివర;

  తిరిగి;
}
