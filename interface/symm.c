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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_نشده COMPLEX
#تعریف_کن SMP_آستانه_MIN 65536.
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DSYMM "
#وإلا
#تعریف_کن ERROR_NAME "SSYMM "
#پایان_اگر
#وإلا
#تعریف_کن SMP_آستانه_MIN 8192.
#اگر_تعریف_نشده GEMM3M
#اگر_تعریف_نشده HEMM
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZSYMM "
#وإلا
#تعریف_کن ERROR_NAME "CSYMM "
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZHEMM "
#وإلا
#تعریف_کن ERROR_NAME "CHEMM "
#پایان_اگر
#پایان_اگر
#وإلا
#اگر_تعریف_نشده HEMM
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZSYMM3M "
#وإلا
#تعریف_کن ERROR_NAME "CSYMM3M "
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZHEMM3M "
#وإلا
#تعریف_کن ERROR_NAME "CHEMM3M "
#پایان_اگر
#پایان_اگر
#پایان_اگر
#پایان_اگر

#اگر_تعریف_نشده GEMM_MULTITHREAD_آستانه
#تعریف_کن GEMM_MULTITHREAD_آستانه 4
#پایان_اگر

#اگر_تعریف_شده DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#پایان_اگر                                                               

#اگر_تعریف_شده SMP
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن نما_2	(BLAS_XDOUBLE | BLAS_REAL)
#وإلاإذا defined(DOUBLE)
#تعریف_کن نما_2	(BLAS_DOUBLE  | BLAS_REAL)
#وإلا
#تعریف_کن نما_2	(BLAS_SINGLE  | BLAS_REAL)
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن نما_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#تعریف_کن نما_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#وإلا
#تعریف_کن نما_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#پایان_اگر
#پایان_اگر
#پایان_اگر

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_تعریف_نشده GEMM3M
#اگر_تعریف_نشده HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#پایان_اگر
#وإلا
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#پایان_اگر
#پایان_اگر
#وإلا
#اگر_تعریف_نشده HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#پایان_اگر
#وإلا
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#پایان_اگر
#پایان_اگر
#پایان_اگر
};

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t آرگومانها;

  FLOAT *میانگیر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#پایان_اگر
# إذا defined(SMP)
  مزدوج MN;
#پایان_اگر
  blasint info;
  صحيح side;
  صحيح uplo;

  PRINT_DEBUG_NAME;

  آرگومانها.alpha = (فراغ *)alpha;
  آرگومانها.beta  = (فراغ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  آرگومانها.m = *M;
  آرگومانها.n = *N;

  آرگومانها.c = (فراغ *)c;
  آرگومانها.ldc = *ldC;

  info = 0;

  إذا (آرگومانها.ldc < MAX(1, آرگومانها.m)) info = 12;

  إذا (!side) {
    آرگومانها.a = (فراغ *)a;
    آرگومانها.b = (فراغ *)b;

    آرگومانها.lda = *ldA;
    آرگومانها.ldb = *ldB;

    إذا (آرگومانها.ldb < MAX(1, آرگومانها.m)) info =  9;
    إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  7;

  } وإلا {
    آرگومانها.a = (فراغ *)b;
    آرگومانها.b = (فراغ *)a;

    آرگومانها.lda = *ldB;
    آرگومانها.ldb = *ldA;

  إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  9;
  إذا (آرگومانها.ldb < MAX(1, آرگومانها.n)) info =  7;
  }

  إذا (آرگومانها.n   < 0)              info =  4;
  إذا (آرگومانها.m   < 0)              info =  3;
  إذا (uplo     < 0)              info =  2;
  إذا (side     < 0)              info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#اگر_تعریف_نشده COMPLEX
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
#پایان_اگر

  blas_arg_t آرگومانها;
  صحيح side, uplo;
  blasint info;

  FLOAT *میانگیر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#پایان_اگر
#إذا defined(SMP)
  مزدوج MN;
#پایان_اگر

  PRINT_DEBUG_CNAME;

#اگر_تعریف_نشده COMPLEX
  آرگومانها.alpha = (فراغ *)&alpha;
  آرگومانها.beta  = (فراغ *)&beta;
#وإلا
  آرگومانها.alpha = (فراغ *)alpha;
  آرگومانها.beta  = (فراغ *)beta;
#پایان_اگر

  آرگومانها.c = (فراغ *)c;
  آرگومانها.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Side == Cblas_چپ)  side = 0;
    إذا (Side == Cblas_راست) side = 1;

    إذا (Uplo == Cblas_فوقانی) uplo  = 0;
    إذا (Uplo == Cblas_پایینی) uplo  = 1;

    info = -1;

    آرگومانها.m = m;
    آرگومانها.n = n;

    إذا (آرگومانها.ldc < MAX(1, آرگومانها.m)) info = 12;

    إذا (!side) {
      آرگومانها.a = (فراغ *)a;
      آرگومانها.b = (فراغ *)b;

      آرگومانها.lda = lda;
      آرگومانها.ldb = ldb;

      إذا (آرگومانها.ldb < MAX(1, آرگومانها.m)) info =  9;
      إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  7;

    } وإلا {
      آرگومانها.a = (فراغ *)b;
      آرگومانها.b = (فراغ *)a;

      آرگومانها.lda = ldb;
      آرگومانها.ldb = lda;

      إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  9;
      إذا (آرگومانها.ldb < MAX(1, آرگومانها.n)) info =  7;
    }

    إذا (آرگومانها.n   < 0)              info =  4;
    إذا (آرگومانها.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (ترتیب == Cblas_ردیف_Major) {
    إذا (Side == Cblas_چپ)  side = 1;
    إذا (Side == Cblas_راست) side = 0;

    إذا (Uplo == Cblas_فوقانی) uplo  = 1;
    إذا (Uplo == Cblas_پایینی) uplo  = 0;

    info = -1;

    آرگومانها.m = n;
    آرگومانها.n = m;

    إذا (آرگومانها.ldc < MAX(1, آرگومانها.m)) info = 12;

    إذا (!side) {
      آرگومانها.a = (فراغ *)a;
      آرگومانها.b = (فراغ *)b;

      آرگومانها.lda = lda;
      آرگومانها.ldb = ldb;

      إذا (آرگومانها.ldb < MAX(1, آرگومانها.m)) info =  9;
      إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  7;

    } وإلا {
      آرگومانها.a = (فراغ *)b;
      آرگومانها.b = (فراغ *)a;

      آرگومانها.lda = ldb;
      آرگومانها.ldb = lda;

      إذا (آرگومانها.lda < MAX(1, آرگومانها.m)) info =  9;
      إذا (آرگومانها.ldb < MAX(1, آرگومانها.n)) info =  7;
    }

    إذا (آرگومانها.n   < 0)              info =  4;
    إذا (آرگومانها.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

   إذا (آرگومانها.m == 0 || آرگومانها.n == 0) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#پایان_اگر
)
#پایان_اگر
   إذا (ترتیب == Cblas_ردیف_Major && m == lda && n == ldb && n == ldc)
   {
     إذا (Side == Cblas_چپ && Uplo == Cblas_فوقانی) {
       SSYMM_مستقیم_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
     وإلا إذا (Side == Cblas_چپ && Uplo == Cblas_پایینی) {
       SSYMM_مستقیم_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
   }
#پایان_اگر
#پایان_اگر

#پایان_اگر

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  میانگیر = (FLOAT *)blas_حافظه_alloc(0);
  إذا (!میانگیر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);

#اگر_تعریف_شده SMP
  آرگومانها.common = عدم;
  MN = 2.* (مزدوج) آرگومانها.m * (مزدوج)آرگومانها.m * (مزدوج) آرگومانها.n;
  إذا (MN <= (SMP_آستانه_MIN * (مزدوج) GEMM_MULTITHREAD_آستانه) ) {
	  آرگومانها.nthreads = 1;
  } وإلا {
	 آرگومانها.nthreads = رقم_cpu_avail(3);
  }
  إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

    (symm[(side << 1) | uplo ])(&آرگومانها, عدم, عدم, sa, sb, 0);

#اگر_تعریف_شده SMP
  } وإلا {

#اگر_تعریف_نشده NO_AFFINITY
    nodes = get_رقم_nodes();

    إذا (nodes > 1) {

      آرگومانها.nthreads /= nodes;

      gemm_thread_mn(نما_2, &آرگومانها, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } وإلا {
#پایان_اگر

#اگر_تعریف_نشده USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&آرگومانها, عدم, عدم, sa, sb, 0);

#وإلا

      GEMM_THREAD(نما_2, &آرگومانها, عدم, عدم, symm[(side << 1) | uplo ], sa, sb, آرگومانها.nthreads);

#پایان_اگر

#اگر_تعریف_نشده NO_AFFINITY
    }
#پایان_اگر

  }
#پایان_اگر

 blas_حافظه_free(میانگیر);

  FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE,
		       (!side)? آرگومانها.m * (آرگومانها.m / 2 + آرگومانها.n) : آرگومانها.n * (آرگومانها.m + آرگومانها.n / 2),
		       (!side)? 2 * آرگومانها.m * آرگومانها.m * آرگومانها.n : 2 * آرگومانها.m * آرگومانها.n * آرگومانها.n);

  IDEBUG_انتها;

  أرجع;
}
