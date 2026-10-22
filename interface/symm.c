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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#که_نه_تعریف_شوی COMPLEX
#تعریف_کړه SMP_پوله_MIN 65536.
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DSYMM "
#وإلا
#تعریف_کړه ERROR_NAME "SSYMM "
#د_که_پای
#وإلا
#تعریف_کړه SMP_پوله_MIN 8192.
#که_نه_تعریف_شوی GEMM3M
#که_نه_تعریف_شوی HEMM
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XSYMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZSYMM "
#وإلا
#تعریف_کړه ERROR_NAME "CSYMM "
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XHEMM "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZHEMM "
#وإلا
#تعریف_کړه ERROR_NAME "CHEMM "
#د_که_پای
#د_که_پای
#وإلا
#که_نه_تعریف_شوی HEMM
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XSYMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZSYMM3M "
#وإلا
#تعریف_کړه ERROR_NAME "CSYMM3M "
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XHEMM3M "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZHEMM3M "
#وإلا
#تعریف_کړه ERROR_NAME "CHEMM3M "
#د_که_پای
#د_که_پای
#د_که_پای
#د_که_پای

#که_نه_تعریف_شوی GEMM_MULTITHREAD_پوله
#تعریف_کړه GEMM_MULTITHREAD_پوله 4
#د_که_پای

#که_تعریف_شوی DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);                                
#د_که_پای                                                               

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ډېرپېښ_2	(BLAS_XDOUBLE | BLAS_REAL)
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ډېرپېښ_2	(BLAS_DOUBLE  | BLAS_REAL)
#وإلا
#تعریف_کړه ډېرپېښ_2	(BLAS_SINGLE  | BLAS_REAL)
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ډېرپېښ_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ډېرپېښ_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#وإلا
#تعریف_کړه ډېرپېښ_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#د_که_پای
#د_که_پای
#د_که_پای

ساكن صحيح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#که_نه_تعریف_شوی GEMM3M
#که_نه_تعریف_شوی HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#د_که_پای
#وإلا
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#د_که_پای
#د_که_پای
#وإلا
#که_نه_تعریف_شوی HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#د_که_پای
#وإلا
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#د_که_پای
#د_که_پای
#د_که_پای
};

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t دلیلونه;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#د_که_پای
# إذا defined(SMP)
  مزدوج MN;
#د_که_پای
  blasint info;
  صحيح side;
  صحيح uplo;

  PRINT_DEBUG_NAME;

  دلیلونه.alpha = (فراغ *)alpha;
  دلیلونه.beta  = (فراغ *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  إذا (side_arg  == 'L') side  = 0;
  إذا (side_arg  == 'R') side  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  دلیلونه.m = *M;
  دلیلونه.n = *N;

  دلیلونه.c = (فراغ *)c;
  دلیلونه.ldc = *ldC;

  info = 0;

  إذا (دلیلونه.ldc < MAX(1, دلیلونه.m)) info = 12;

  إذا (!side) {
    دلیلونه.a = (فراغ *)a;
    دلیلونه.b = (فراغ *)b;

    دلیلونه.lda = *ldA;
    دلیلونه.ldb = *ldB;

    إذا (دلیلونه.ldb < MAX(1, دلیلونه.m)) info =  9;
    إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  7;

  } وإلا {
    دلیلونه.a = (فراغ *)b;
    دلیلونه.b = (فراغ *)a;

    دلیلونه.lda = *ldB;
    دلیلونه.ldb = *ldA;

  إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  9;
  إذا (دلیلونه.ldb < MAX(1, دلیلونه.n)) info =  7;
  }

  إذا (دلیلونه.n   < 0)              info =  4;
  إذا (دلیلونه.m   < 0)              info =  3;
  إذا (uplo     < 0)              info =  2;
  إذا (side     < 0)              info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_SIDE Side, تعداد CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#که_نه_تعریف_شوی COMPLEX
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
#د_که_پای

  blas_arg_t دلیلونه;
  صحيح side, uplo;
  blasint info;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#إذا defined(SMP) && !defined(NO_AFFINITY)
  صحيح nodes;
#د_که_پای
#إذا defined(SMP)
  مزدوج MN;
#د_که_پای

  PRINT_DEBUG_CNAME;

#که_نه_تعریف_شوی COMPLEX
  دلیلونه.alpha = (فراغ *)&alpha;
  دلیلونه.beta  = (فراغ *)&beta;
#وإلا
  دلیلونه.alpha = (فراغ *)alpha;
  دلیلونه.beta  = (فراغ *)beta;
#د_که_پای

  دلیلونه.c = (فراغ *)c;
  دلیلونه.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_ستون_Major) {
    إذا (Side == Cblas_کیڼ)  side = 0;
    إذا (Side == Cblas_ښی) side = 1;

    إذا (Uplo == Cblas_پاسنی) uplo  = 0;
    إذا (Uplo == Cblas_ښکتنی) uplo  = 1;

    info = -1;

    دلیلونه.m = m;
    دلیلونه.n = n;

    إذا (دلیلونه.ldc < MAX(1, دلیلونه.m)) info = 12;

    إذا (!side) {
      دلیلونه.a = (فراغ *)a;
      دلیلونه.b = (فراغ *)b;

      دلیلونه.lda = lda;
      دلیلونه.ldb = ldb;

      إذا (دلیلونه.ldb < MAX(1, دلیلونه.m)) info =  9;
      إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  7;

    } وإلا {
      دلیلونه.a = (فراغ *)b;
      دلیلونه.b = (فراغ *)a;

      دلیلونه.lda = ldb;
      دلیلونه.ldb = lda;

      إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  9;
      إذا (دلیلونه.ldb < MAX(1, دلیلونه.n)) info =  7;
    }

    إذا (دلیلونه.n   < 0)              info =  4;
    إذا (دلیلونه.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (ترتیب == Cblas_لیکه_Major) {
    إذا (Side == Cblas_کیڼ)  side = 1;
    إذا (Side == Cblas_ښی) side = 0;

    إذا (Uplo == Cblas_پاسنی) uplo  = 1;
    إذا (Uplo == Cblas_ښکتنی) uplo  = 0;

    info = -1;

    دلیلونه.m = n;
    دلیلونه.n = m;

    إذا (دلیلونه.ldc < MAX(1, دلیلونه.m)) info = 12;

    إذا (!side) {
      دلیلونه.a = (فراغ *)a;
      دلیلونه.b = (فراغ *)b;

      دلیلونه.lda = lda;
      دلیلونه.ldb = ldb;

      إذا (دلیلونه.ldb < MAX(1, دلیلونه.m)) info =  9;
      إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  7;

    } وإلا {
      دلیلونه.a = (فراغ *)b;
      دلیلونه.b = (فراغ *)a;

      دلیلونه.lda = ldb;
      دلیلونه.ldb = lda;

      إذا (دلیلونه.lda < MAX(1, دلیلونه.m)) info =  9;
      إذا (دلیلونه.ldb < MAX(1, دلیلونه.n)) info =  7;
    }

    إذا (دلیلونه.n   < 0)              info =  4;
    إذا (دلیلونه.m   < 0)              info =  3;
    إذا (uplo     < 0)              info =  2;
    إذا (side     < 0)              info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

   إذا (دلیلونه.m == 0 || دلیلونه.n == 0) أرجع;
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_مستقیم)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#د_که_پای
)
#د_که_پای
   إذا (ترتیب == Cblas_لیکه_Major && m == lda && n == ldb && n == ldc)
   {
     إذا (Side == Cblas_کیڼ && Uplo == Cblas_پاسنی) {
       SSYMM_مستقیم_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
     وإلا إذا (Side == Cblas_کیڼ && Uplo == Cblas_ښکتنی) {
       SSYMM_مستقیم_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); أرجع;
     }
   }
#د_که_پای
#د_که_پای

#د_که_پای

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  بفر = (FLOAT *)blas_حافظه_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_دننه_واټن_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_دننه_واټن_B);

#که_تعریف_شوی SMP
  دلیلونه.common = عدم;
  MN = 2.* (مزدوج) دلیلونه.m * (مزدوج)دلیلونه.m * (مزدوج) دلیلونه.n;
  إذا (MN <= (SMP_پوله_MIN * (مزدوج) GEMM_MULTITHREAD_پوله) ) {
	  دلیلونه.nthreads = 1;
  } وإلا {
	 دلیلونه.nthreads = عدد_cpu_avail(3);
  }
  إذا (دلیلونه.nthreads == 1) {
#د_که_پای

    (symm[(side << 1) | uplo ])(&دلیلونه, عدم, عدم, sa, sb, 0);

#که_تعریف_شوی SMP
  } وإلا {

#که_نه_تعریف_شوی NO_AFFINITY
    nodes = get_عدد_nodes();

    إذا (nodes > 1) {

      دلیلونه.nthreads /= nodes;

      gemm_thread_mn(ډېرپېښ_2, &دلیلونه, عدم, عدم,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } وإلا {
#د_که_پای

#که_نه_تعریف_شوی USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&دلیلونه, عدم, عدم, sa, sb, 0);

#وإلا

      GEMM_THREAD(ډېرپېښ_2, &دلیلونه, عدم, عدم, symm[(side << 1) | uplo ], sa, sb, دلیلونه.nthreads);

#د_که_پای

#که_نه_تعریف_شوی NO_AFFINITY
    }
#د_که_پای

  }
#د_که_پای

 blas_حافظه_free(بفر);

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE * COMPSIZE,
		       (!side)? دلیلونه.m * (دلیلونه.m / 2 + دلیلونه.n) : دلیلونه.n * (دلیلونه.m + دلیلونه.n / 2),
		       (!side)? 2 * دلیلونه.m * دلیلونه.m * دلیلونه.n : 2 * دلیلونه.m * دلیلونه.n * دلیلونه.n);

  IDEBUG_پای_ټکی;

  أرجع;
}
