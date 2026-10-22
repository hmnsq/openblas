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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_بنا_تعریف COMPLEX
#تعریف_کرو SMP_حد_MIN 65536.
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DSYMM "
#وإلا
#تعریف_کرو ERROR_NAME "SSYMM "
#جے_انت
#وإلا
#تعریف_کرو SMP_حد_MIN 8192.
#جے_بنا_تعریف GEMM3M
#جے_بنا_تعریف HEMM
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYMM "
#وإلا
#تعریف_کرو ERROR_NAME "CSYMM "
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHEMM "
#وإلا
#تعریف_کرو ERROR_NAME "CHEMM "
#جے_انت
#جے_انت
#وإلا
#جے_بنا_تعریف HEMM
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYMM3M "
#وإلا
#تعریف_کرو ERROR_NAME "CSYMM3M "
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHEMM3M "
#وإلا
#تعریف_کرو ERROR_NAME "CHEMM3M "
#جے_انت
#جے_انت
#جے_انت
#جے_انت

#جے_بنا_تعریف GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#جے_انت

#جے_تعریف_شدہ DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#جے_انت                                                               

#جے_تعریف_شدہ SMP
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو منوال_2	(BLAS_XDOUBLE | BLAS_REAL)
#وإلاإذا defined(DOUBLE)
#تعریف_کرو منوال_2	(BLAS_DOUBLE  | BLAS_REAL)
#وإلا
#تعریف_کرو منوال_2	(BLAS_SINGLE  | BLAS_REAL)
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو منوال_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#تعریف_کرو منوال_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#وإلا
#تعریف_کرو منوال_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#جے_انت
#جے_انت
#جے_انت

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جے_بنا_تعریف GEMM3M
#جے_بنا_تعریف HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#جے_انت
#وإلا
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#جے_انت
#جے_انت
#وإلا
#جے_بنا_تعریف HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#جے_انت
#وإلا
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#جے_انت
#جے_انت
#جے_انت
};

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t دلیلاں;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#جے_انت
# إذا defined(SMP)
  مزدوج MN;
#جے_انت
  blasint info;
  صحيح side;
  صحيح uplo;

  PRINT_DEBUG_NAME;

  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  دلیلاں.m = *M;
  دلیلاں.n = *N;

  دلیلاں.c = (فراغ *)c;
  دلیلاں.ldc = *ldC;

  info = 0;

  إذا (دلیلاں.ldc < MAX(1, دلیلاں.m)) info = 12;

  إذا (!side) {
    دلیلاں.a = (فراغ *)a;
    دلیلاں.b = (فراغ *)b;

    دلیلاں.lda = *ldA;
    دلیلاں.ldb = *ldB;

    إذا (دلیلاں.ldb < MAX(1, دلیلاں.m)) info =  9;
    إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  7;

  } وإلا {
    دلیلاں.a = (فراغ *)b;
    دلیلاں.b = (فراغ *)a;

    دلیلاں.lda = *ldB;
    دلیلاں.ldb = *ldA;

  إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  9;
  إذا (دلیلاں.ldb < MAX(1, دلیلاں.n)) info =  7;
  }

  إذا (دلیلاں.n   < 0)              info =  4;
  إذا (دلیلاں.m   < 0)              info =  3;
  إذا (uplo     < 0)              info =  2;
  إذا (side     < 0)              info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#جے_بنا_تعریف COMPLEX
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
#جے_انت

  blas_arg_t دلیلاں;
  صحيح side, uplo;
  blasint info;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#جے_انت
#إذا defined(SMP)
  مزدوج MN;
#جے_انت

  PRINT_DEBUG_CNAME;

#جے_بنا_تعریف COMPLEX
  دلیلاں.alpha = (فراغ *)&alpha;
  دلیلاں.beta  = (فراغ *)&beta;
#وإلا
  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;
#جے_انت

  دلیلاں.c = (فراغ *)c;
  دلیلاں.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {
    إذا (Side == Cblas_بایاں)  side = 0;
    إذا (Side == Cblas_دایاں) side = 1;

    إذا (Uplo == Cblas_اتلا) uplo  = 0;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 1;

    info = -1;

    دلیلاں.m = m;
    دلیلاں.n = n;

    إذا (دلیلاں.ldc < MAX(1, دلیلاں.m)) info = 12;

    إذا (!side) {
      دلیلاں.a = (فراغ *)a;
      دلیلاں.b = (فراغ *)b;

      دلیلاں.lda = lda;
      دلیلاں.ldb = ldb;

      إذا (دلیلاں.ldb < MAX(1, دلیلاں.m)) info =  9;
      إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  7;

    } وإلا {
      دلیلاں.a = (فراغ *)b;
      دلیلاں.b = (فراغ *)a;

      دلیلاں.lda = ldb;
      دلیلاں.ldb = lda;

      إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  9;
      إذا (دلیلاں.ldb < MAX(1, دلیلاں.n)) info =  7;
    }

    إذا (دلیلاں.n   < 0)              info =  4;
    إذا (دلیلاں.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (ترتیب == Cblas_صف_Major) {
    إذا (Side == Cblas_بایاں)  side = 1;
    إذا (Side == Cblas_دایاں) side = 0;

    إذا (Uplo == Cblas_اتلا) uplo  = 1;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 0;

    info = -1;

    دلیلاں.m = n;
    دلیلاں.n = m;

    إذا (دلیلاں.ldc < MAX(1, دلیلاں.m)) info = 12;

    إذا (!side) {
      دلیلاں.a = (فراغ *)a;
      دلیلاں.b = (فراغ *)b;

      دلیلاں.lda = lda;
      دلیلاں.ldb = ldb;

      إذا (دلیلاں.ldb < MAX(1, دلیلاں.m)) info =  9;
      إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  7;

    } وإلا {
      دلیلاں.a = (فراغ *)b;
      دلیلاں.b = (فراغ *)a;

      دلیلاں.lda = ldb;
      دلیلاں.ldb = lda;

      إذا (دلیلاں.lda < MAX(1, دلیلاں.m)) info =  9;
      إذا (دلیلاں.ldb < MAX(1, دلیلاں.n)) info =  7;
    }

    إذا (دلیلاں.n   < 0)              info =  4;
    إذا (دلیلاں.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

   إذا (دلیلاں.m == 0 || دلیلاں.n == 0) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جے_انت
)
#جے_انت
   إذا (ترتیب == Cblas_صف_Major && m == lda && n == ldb && n == ldc)
   {
     إذا (Side == Cblas_بایاں && Uplo == Cblas_اتلا) {
       SSYMM_سدھا_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
     وإلا إذا (Side == Cblas_بایاں && Uplo == Cblas_ہیٹھلا) {
       SSYMM_سدھا_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
   }
#جے_انت
#جے_انت

#جے_انت

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وتھ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ناپ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وتھ_B);

#جے_تعریف_شدہ SMP
  دلیلاں.common = عدم;
  MN = 2.* (مزدوج) دلیلاں.m * (مزدوج)دلیلاں.m * (مزدوج) دلیلاں.n;
  إذا (MN <= (SMP_حد_MIN * (مزدوج) GEMM_MULTITHREAD_حد) ) {
	  دلیلاں.nthreads = 1;
  } وإلا {
	 دلیلاں.nthreads = عدد_cpu_avail(3);
  }
  إذا (دلیلاں.nthreads == 1) {
#جے_انت

    (symm[(side << 1) | uplo ])(&دلیلاں, عدم, عدم, sa, sb, 0);

#جے_تعریف_شدہ SMP
  } وإلا {

#جے_بنا_تعریف NO_AFFINITY
    nodes = get_عدد_nodes();

    إذا (nodes > 1) {

      دلیلاں.nthreads /= nodes;

      gemm_thread_mn(منوال_2, &دلیلاں, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } وإلا {
#جے_انت

#جے_بنا_تعریف USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&دلیلاں, عدم, عدم, sa, sb, 0);

#وإلا

      GEMM_THREAD(منوال_2, &دلیلاں, عدم, عدم, symm[(side << 1) | uplo ], sa, sb, دلیلاں.nthreads);

#جے_انت

#جے_بنا_تعریف NO_AFFINITY
    }
#جے_انت

  }
#جے_انت

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_سرا(COMPSIZE * COMPSIZE,
		       (!side)? دلیلاں.m * (دلیلاں.m / 2 + دلیلاں.n) : دلیلاں.n * (دلیلاں.m + دلیلاں.n / 2),
		       (!side)? 2 * دلیلاں.m * دلیلاں.m * دلیلاں.n : 2 * دلیلاں.m * دلیلاں.n * دلیلاں.n);

  IDEBUG_سرا;

  أرجع;
}
