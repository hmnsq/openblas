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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#وصف_ڪريو SMP_حد_MIN 65536.
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QSYMM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DSYMM "
#وإلا
#وصف_ڪريو ERROR_NAME "SSYMM "
#جيڪڏهن_پورو
#وإلا
#وصف_ڪريو SMP_حد_MIN 8192.
#جيڪڏهن_وصف_نه_ٿيل GEMM3M
#جيڪڏهن_وصف_نه_ٿيل HEMM
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZSYMM "
#وإلا
#وصف_ڪريو ERROR_NAME "CSYMM "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHEMM "
#وإلا
#وصف_ڪريو ERROR_NAME "CHEMM "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل HEMM
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZSYMM3M "
#وإلا
#وصف_ڪريو ERROR_NAME "CSYMM3M "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHEMM3M "
#وإلا
#وصف_ڪريو ERROR_NAME "CHEMM3M "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل GEMM_MULTITHREAD_حد
#وصف_ڪريو GEMM_MULTITHREAD_حد 4
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#جيڪڏهن_پورو                                                               

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو گهڻائي_2	(BLAS_XDOUBLE | BLAS_REAL)
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو گهڻائي_2	(BLAS_DOUBLE  | BLAS_REAL)
#وإلا
#وصف_ڪريو گهڻائي_2	(BLAS_SINGLE  | BLAS_REAL)
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو گهڻائي_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو گهڻائي_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#وإلا
#وصف_ڪريو گهڻائي_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جيڪڏهن_وصف_نه_ٿيل GEMM3M
#جيڪڏهن_وصف_نه_ٿيل HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#جيڪڏهن_پورو
#وإلا
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_نه_ٿيل HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#جيڪڏهن_پورو
#وإلا
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو
};

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t دليلون;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#جيڪڏهن_پورو
# إذا defined(SMP)
  مزدوج MN;
#جيڪڏهن_پورو
  blasint info;
  صحيح side;
  صحيح uplo;

  PRINT_DEBUG_NAME;

  دليلون.alpha = (فراغ *)alpha;
  دليلون.beta  = (فراغ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  دليلون.m = *M;
  دليلون.n = *N;

  دليلون.c = (فراغ *)c;
  دليلون.ldc = *ldC;

  info = 0;

  إذا (دليلون.ldc < MAX(1, دليلون.m)) info = 12;

  إذا (!side) {
    دليلون.a = (فراغ *)a;
    دليلون.b = (فراغ *)b;

    دليلون.lda = *ldA;
    دليلون.ldb = *ldB;

    إذا (دليلون.ldb < MAX(1, دليلون.m)) info =  9;
    إذا (دليلون.lda < MAX(1, دليلون.m)) info =  7;

  } وإلا {
    دليلون.a = (فراغ *)b;
    دليلون.b = (فراغ *)a;

    دليلون.lda = *ldB;
    دليلون.ldb = *ldA;

  إذا (دليلون.lda < MAX(1, دليلون.m)) info =  9;
  إذا (دليلون.ldb < MAX(1, دليلون.n)) info =  7;
  }

  إذا (دليلون.n   < 0)              info =  4;
  إذا (دليلون.m   < 0)              info =  3;
  إذا (uplo     < 0)              info =  2;
  إذا (side     < 0)              info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#وإلا
	   فراغ *valpha,
	   فراغ *va, blasint lda,
	   فراغ *vb, blasint ldb,
	   فراغ *vbeta,
	   فراغ *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#جيڪڏهن_پورو

  blas_arg_t دليلون;
  صحيح side, uplo;
  blasint info;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#جيڪڏهن_پورو
#إذا defined(SMP)
  مزدوج MN;
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  دليلون.alpha = (فراغ *)&alpha;
  دليلون.beta  = (فراغ *)&beta;
#وإلا
  دليلون.alpha = (فراغ *)alpha;
  دليلون.beta  = (فراغ *)beta;
#جيڪڏهن_پورو

  دليلون.c = (فراغ *)c;
  دليلون.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (Side == Cblas_کاٻو)  side = 0;
    إذا (Side == Cblas_ساڄو) side = 1;

    إذا (Uplo == Cblas_مٿانهون) uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 1;

    info = -1;

    دليلون.m = m;
    دليلون.n = n;

    إذا (دليلون.ldc < MAX(1, دليلون.m)) info = 12;

    إذا (!side) {
      دليلون.a = (فراغ *)a;
      دليلون.b = (فراغ *)b;

      دليلون.lda = lda;
      دليلون.ldb = ldb;

      إذا (دليلون.ldb < MAX(1, دليلون.m)) info =  9;
      إذا (دليلون.lda < MAX(1, دليلون.m)) info =  7;

    } وإلا {
      دليلون.a = (فراغ *)b;
      دليلون.b = (فراغ *)a;

      دليلون.lda = ldb;
      دليلون.ldb = lda;

      إذا (دليلون.lda < MAX(1, دليلون.m)) info =  9;
      إذا (دليلون.ldb < MAX(1, دليلون.n)) info =  7;
    }

    إذا (دليلون.n   < 0)              info =  4;
    إذا (دليلون.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (Side == Cblas_کاٻو)  side = 1;
    إذا (Side == Cblas_ساڄو) side = 0;

    إذا (Uplo == Cblas_مٿانهون) uplo  = 1;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 0;

    info = -1;

    دليلون.m = n;
    دليلون.n = m;

    إذا (دليلون.ldc < MAX(1, دليلون.m)) info = 12;

    إذا (!side) {
      دليلون.a = (فراغ *)a;
      دليلون.b = (فراغ *)b;

      دليلون.lda = lda;
      دليلون.ldb = ldb;

      إذا (دليلون.ldb < MAX(1, دليلون.m)) info =  9;
      إذا (دليلون.lda < MAX(1, دليلون.m)) info =  7;

    } وإلا {
      دليلون.a = (فراغ *)b;
      دليلون.b = (فراغ *)a;

      دليلون.lda = ldb;
      دليلون.ldb = lda;

      إذا (دليلون.lda < MAX(1, دليلون.m)) info =  9;
      إذا (دليلون.ldb < MAX(1, دليلون.n)) info =  7;
    }

    إذا (دليلون.n   < 0)              info =  4;
    إذا (دليلون.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

   إذا (دليلون.m == 0 || دليلون.n == 0) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_سڌو)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جيڪڏهن_پورو
)
#جيڪڏهن_پورو
   إذا (ترتيب == Cblas_صف_Major && m == lda && n == ldb && n == ldc)
   {
     إذا (Side == Cblas_کاٻو && Uplo == Cblas_مٿانهون) {
       SSYMM_سڌو_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
     وإلا إذا (Side == Cblas_کاٻو && Uplo == Cblas_هيٺاهون) {
       SSYMM_سڌو_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
   }
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_پورو

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  بفر = (FLOAT *)blas_يادگيري_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);

#جيڪڏهن_وصف_ٿيل SMP
  دليلون.common = عدم;
  MN = 2.* (مزدوج) دليلون.m * (مزدوج)دليلون.m * (مزدوج) دليلون.n;
  إذا (MN <= (SMP_حد_MIN * (مزدوج) GEMM_MULTITHREAD_حد) ) {
	  دليلون.nthreads = 1;
  } وإلا {
	 دليلون.nthreads = انگ_cpu_avail(3);
  }
  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    (symm[(side << 1) | uplo ])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

#جيڪڏهن_وصف_نه_ٿيل NO_AFFINITY
    nodes = get_انگ_nodes();

    إذا (nodes > 1) {

      دليلون.nthreads /= nodes;

      gemm_thread_mn(گهڻائي_2, &دليلون, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } وإلا {
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&دليلون, عدم, عدم, sa, sb, 0);

#وإلا

      GEMM_THREAD(گهڻائي_2, &دليلون, عدم, عدم, symm[(side << 1) | uplo ], sa, sb, دليلون.nthreads);

#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل NO_AFFINITY
    }
#جيڪڏهن_پورو

  }
#جيڪڏهن_پورو

 blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE,
		       (!side)? دليلون.m * (دليلون.m / 2 + دليلون.n) : دليلون.n * (دليلون.m + دليلون.n / 2),
		       (!side)? 2 * دليلون.m * دليلون.m * دليلون.n : 2 * دليلون.m * دليلون.n * دليلون.n);

  IDEBUG_پڇاڙي;

  أرجع;
}
