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

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو SMP_حد_MIN 65536.
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QSYMM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DSYMM "
#ورنہ
#تعریف_کرو ERROR_NAME "SSYMM "
#اگر_ختم
#ورنہ
#تعریف_کرو SMP_حد_MIN 8192.
#اگر_غیر_تعریف_شدہ GEMM3M
#اگر_غیر_تعریف_شدہ HEMM
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYMM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYMM "
#ورنہ
#تعریف_کرو ERROR_NAME "CSYMM "
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHEMM "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHEMM "
#ورنہ
#تعریف_کرو ERROR_NAME "CHEMM "
#اگر_ختم
#اگر_ختم
#ورنہ
#اگر_غیر_تعریف_شدہ HEMM
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYMM3M "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYMM3M "
#ورنہ
#تعریف_کرو ERROR_NAME "CSYMM3M "
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHEMM3M "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHEMM3M "
#ورنہ
#تعریف_کرو ERROR_NAME "CHEMM3M "
#اگر_ختم
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر_غیر_تعریف_شدہ GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#اگر_ختم

#اگر_تعریف_شدہ DYNAMIC_ARCH
خارجی حرف* gotoblas_corename(خالی);                                
#اگر_ختم                                                               

#اگر_تعریف_شدہ SMP
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو منوال_2	(BLAS_XDOUBLE | BLAS_REAL)
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو منوال_2	(BLAS_DOUBLE  | BLAS_REAL)
#ورنہ
#تعریف_کرو منوال_2	(BLAS_SINGLE  | BLAS_REAL)
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو منوال_2	(BLAS_XDOUBLE | BLAS_COMPLEX)
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو منوال_2	(BLAS_DOUBLE  | BLAS_COMPLEX)
#ورنہ
#تعریف_کرو منوال_2	(BLAS_SINGLE  | BLAS_COMPLEX)
#اگر_ختم
#اگر_ختم
#اگر_ختم

ساکن صحیح (*symm[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_غیر_تعریف_شدہ GEMM3M
#اگر_غیر_تعریف_شدہ HEMM
  SYMM_LU, SYMM_LL, SYMM_RU, SYMM_RL,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM_THREAD_LU, SYMM_THREAD_LL, SYMM_THREAD_RU, SYMM_THREAD_RL,
#اگر_ختم
#ورنہ
  HEMM_LU, HEMM_LL, HEMM_RU, HEMM_RL,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM_THREAD_LU, HEMM_THREAD_LL, HEMM_THREAD_RU, HEMM_THREAD_RL,
#اگر_ختم
#اگر_ختم
#ورنہ
#اگر_غیر_تعریف_شدہ HEMM
  SYMM3M_LU, SYMM3M_LL, SYMM3M_RU, SYMM3M_RL,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYMM3M_THREAD_LU, SYMM3M_THREAD_LL, SYMM3M_THREAD_RU, SYMM3M_THREAD_RL,
#اگر_ختم
#ورنہ
  HEMM3M_LU, HEMM3M_LL, HEMM3M_RU, HEMM3M_RL,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HEMM3M_THREAD_LU, HEMM3M_THREAD_LL, HEMM3M_THREAD_RU, HEMM3M_THREAD_RL,
#اگر_ختم
#اگر_ختم
#اگر_ختم
};

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *SIDE, حرف *UPLO,
         blasint *M, blasint *N,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *b, blasint *ldB,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف side_arg  = *SIDE;
  حرف uplo_arg  = *UPLO;

  blas_arg_t دلائل;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#اگر defined(SMP) && !defined(NO_AFFINITY)
  صحیح nodes;
#اگر_ختم
# اگر defined(SMP)
  دوہرا MN;
#اگر_ختم
  blasint info;
  صحیح side;
  صحیح uplo;

  PRINT_DEBUG_NAME;

  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;

  TOUPPER(side_arg);
  TOUPPER(uplo_arg);

  side  = -1;
  uplo  = -1;

  اگر (side_arg  == 'L') side  = 0;
  اگر (side_arg  == 'R') side  = 1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  دلائل.m = *M;
  دلائل.n = *N;

  دلائل.c = (خالی *)c;
  دلائل.ldc = *ldC;

  info = 0;

  اگر (دلائل.ldc < MAX(1, دلائل.m)) info = 12;

  اگر (!side) {
    دلائل.a = (خالی *)a;
    دلائل.b = (خالی *)b;

    دلائل.lda = *ldA;
    دلائل.ldb = *ldB;

    اگر (دلائل.ldb < MAX(1, دلائل.m)) info =  9;
    اگر (دلائل.lda < MAX(1, دلائل.m)) info =  7;

  } ورنہ {
    دلائل.a = (خالی *)b;
    دلائل.b = (خالی *)a;

    دلائل.lda = *ldB;
    دلائل.ldb = *ldA;

  اگر (دلائل.lda < MAX(1, دلائل.m)) info =  9;
  اگر (دلائل.ldb < MAX(1, دلائل.n)) info =  7;
  }

  اگر (دلائل.n   < 0)              info =  4;
  اگر (دلائل.m   < 0)              info =  3;
  اگر (uplo     < 0)              info =  2;
  اگر (side     < 0)              info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_SIDE Side, شمار CBLAS_UPLO Uplo,
	   blasint m, blasint n,
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT alpha,
	   FLOAT *a, blasint lda,
	   FLOAT *b, blasint ldb,
	   FLOAT beta,
	   FLOAT *c, blasint ldc) {
#ورنہ
	   خالی *valpha,
	   خالی *va, blasint lda,
	   خالی *vb, blasint ldb,
	   خالی *vbeta,
	   خالی *vc, blasint ldc) {
  FLOAT *alpha = (FLOAT*) valpha;
  FLOAT *beta  = (FLOAT*) vbeta;
  FLOAT *a = (FLOAT*) va;
  FLOAT *b = (FLOAT*) vb;
  FLOAT *c = (FLOAT*) vc;	   
#اگر_ختم

  blas_arg_t دلائل;
  صحیح side, uplo;
  blasint info;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#اگر defined(SMP) && !defined(NO_AFFINITY)
  صحیح nodes;
#اگر_ختم
#اگر defined(SMP)
  دوہرا MN;
#اگر_ختم

  PRINT_DEBUG_CNAME;

#اگر_غیر_تعریف_شدہ COMPLEX
  دلائل.alpha = (خالی *)&alpha;
  دلائل.beta  = (خالی *)&beta;
#ورنہ
  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;
#اگر_ختم

  دلائل.c = (خالی *)c;
  دلائل.ldc = ldc;

  side  = -1;
  uplo  = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    اگر (Side == Cblas_بایاں)  side = 0;
    اگر (Side == Cblas_دایاں) side = 1;

    اگر (Uplo == Cblas_بالائی) uplo  = 0;
    اگر (Uplo == Cblas_زیریں) uplo  = 1;

    info = -1;

    دلائل.m = m;
    دلائل.n = n;

    اگر (دلائل.ldc < MAX(1, دلائل.m)) info = 12;

    اگر (!side) {
      دلائل.a = (خالی *)a;
      دلائل.b = (خالی *)b;

      دلائل.lda = lda;
      دلائل.ldb = ldb;

      اگر (دلائل.ldb < MAX(1, دلائل.m)) info =  9;
      اگر (دلائل.lda < MAX(1, دلائل.m)) info =  7;

    } ورنہ {
      دلائل.a = (خالی *)b;
      دلائل.b = (خالی *)a;

      دلائل.lda = ldb;
      دلائل.ldb = lda;

      اگر (دلائل.lda < MAX(1, دلائل.m)) info =  9;
      اگر (دلائل.ldb < MAX(1, دلائل.n)) info =  7;
    }

    اگر (دلائل.n   < 0)              info =  4;
    اگر (دلائل.m   < 0)              info =  3;
    اگر (uplo     < 0)              info =  2;
    اگر (side     < 0)              info =  1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    اگر (Side == Cblas_بایاں)  side = 1;
    اگر (Side == Cblas_دایاں) side = 0;

    اگر (Uplo == Cblas_بالائی) uplo  = 1;
    اگر (Uplo == Cblas_زیریں) uplo  = 0;

    info = -1;

    دلائل.m = n;
    دلائل.n = m;

    اگر (دلائل.ldc < MAX(1, دلائل.m)) info = 12;

    اگر (!side) {
      دلائل.a = (خالی *)a;
      دلائل.b = (خالی *)b;

      دلائل.lda = lda;
      دلائل.ldb = ldb;

      اگر (دلائل.ldb < MAX(1, دلائل.m)) info =  9;
      اگر (دلائل.lda < MAX(1, دلائل.m)) info =  7;

    } ورنہ {
      دلائل.a = (خالی *)b;
      دلائل.b = (خالی *)a;

      دلائل.lda = ldb;
      دلائل.ldb = lda;

      اگر (دلائل.lda < MAX(1, دلائل.m)) info =  9;
      اگر (دلائل.ldb < MAX(1, دلائل.n)) info =  7;
    }

    اگر (دلائل.n   < 0)              info =  4;
    اگر (دلائل.m   < 0)              info =  3;
    اگر (uplo     < 0)              info =  2;
    اگر (side     < 0)              info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

   اگر (دلائل.m == 0 || دلائل.n == 0) واپس;
#اگر !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_ARM64) && (defined(USE_SSYMM_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#اگر_ختم
)
#اگر_ختم
   اگر (ترتیب == Cblas_صف_Major && m == lda && n == ldb && n == ldc)
   {
     اگر (Side == Cblas_بایاں && Uplo == Cblas_بالائی) {
       SSYMM_براہ_راست_ALPHA_BETA_LU(m, n, alpha, a, lda, b, ldb, beta, c, ldc); واپس;
     }
     ورنہ اگر (Side == Cblas_بایاں && Uplo == Cblas_زیریں) {
       SSYMM_براہ_راست_ALPHA_BETA_LL(m, n, alpha, a, lda, b, ldb, beta, c, ldc); واپس;
     }
   }
#اگر_ختم
#اگر_ختم

#اگر_ختم

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);

#اگر_تعریف_شدہ SMP
  دلائل.common = تہی;
  MN = 2.* (دوہرا) دلائل.m * (دوہرا)دلائل.m * (دوہرا) دلائل.n;
  اگر (MN <= (SMP_حد_MIN * (دوہرا) GEMM_MULTITHREAD_حد) ) {
	  دلائل.nthreads = 1;
  } ورنہ {
	 دلائل.nthreads = عدد_cpu_avail(3);
  }
  اگر (دلائل.nthreads == 1) {
#اگر_ختم

    (symm[(side << 1) | uplo ])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP
  } ورنہ {

#اگر_غیر_تعریف_شدہ NO_AFFINITY
    nodes = get_عدد_nodes();

    اگر (nodes > 1) {

      دلائل.nthreads /= nodes;

      gemm_thread_mn(منوال_2, &دلائل, تہی, تہی,
		     symm[4 | (side << 1) | uplo ], sa, sb, nodes);

    } ورنہ {
#اگر_ختم

#اگر_غیر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3

      (symm[4 | (side << 1) | uplo ])(&دلائل, تہی, تہی, sa, sb, 0);

#ورنہ

      GEMM_THREAD(منوال_2, &دلائل, تہی, تہی, symm[(side << 1) | uplo ], sa, sb, دلائل.nthreads);

#اگر_ختم

#اگر_غیر_تعریف_شدہ NO_AFFINITY
    }
#اگر_ختم

  }
#اگر_ختم

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE,
		       (!side)? دلائل.m * (دلائل.m / 2 + دلائل.n) : دلائل.n * (دلائل.m + دلائل.n / 2),
		       (!side)? 2 * دلائل.m * دلائل.m * دلائل.n : 2 * دلائل.m * دلائل.n * دلائل.n);

  IDEBUG_انتہا;

  واپس;
}
