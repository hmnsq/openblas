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
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QGESV"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DGESV"
#وإلا
#تعریف_کن ERROR_NAME "SGESV"
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XGESV"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZGESV"
#وإلا
#تعریف_کن ERROR_NAME "CGESV"
#پایان_اگر
#پایان_اگر

OPENBLAS_EXPORT
فراغ NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t آرگومانها;

  blasint info;
  FLOAT *میانگیر;
#اگر_تعریف_شده PPC440
  خارجي
#پایان_اگر
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  آرگومانها.m    = *N;
  آرگومانها.n    = *NRHS;
  آرگومانها.a    = (فراغ *)a;
  آرگومانها.lda  = *ldA;
  آرگومانها.b    = (فراغ *)b;
  آرگومانها.ldb  = *ldB;
  آرگومانها.c    = (فراغ *)ipiv;

  info  = 0;
  إذا (آرگومانها.ldb < MAX(1,آرگومانها.m)) info = 7;
  إذا (آرگومانها.lda < MAX(1,آرگومانها.m)) info = 4;
  إذا (آرگومانها.n   < 0)             info = 2;
  إذا (آرگومانها.m   < 0)             info = 1;

  إذا (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    أرجع;
  }

  آرگومانها.alpha = عدم;
  آرگومانها.beta  = عدم;

  *Info = 0;

  إذا (آرگومانها.m == 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#اگر_تعریف_نشده PPC440
  میانگیر = (FLOAT *)blas_حافظه_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);
#پایان_اگر

#اگر_تعریف_شده SMP
  آرگومانها.common = عدم;

#إذا defined(_WIN64) && defined(_M_ARM64)
  #اگر_تعریف_شده COMPLEX
    إذا (آرگومانها.m * آرگومانها.m <= 300) 
  #وإلا
    إذا (آرگومانها.m * آرگومانها.m <= 500) 
  #پایان_اگر
      آرگومانها.nthreads = 1;
  وإلا إذا (آرگومانها.m * آرگومانها.m <= 1000)
    آرگومانها.nthreads = 4;
  وإلا
    آرگومانها.nthreads = رقم_cpu_avail(4);
#وإلا
  #اگر_تعریف_نشده DOUBLE
    إذا (آرگومانها.m * آرگومانها.m < 40000)
  #وإلا
    إذا (آرگومانها.m * آرگومانها.m < 10000)
  #پایان_اگر
      آرگومانها.nthreads = 1;
    وإلا
      آرگومانها.nthreads = رقم_cpu_avail(4);
#پایان_اگر

  إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

    آرگومانها.n    = *N;
    info = GETRF_SINGLE(&آرگومانها, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS >0){
      آرگومانها.n    = *NRHS;
      GETRS_N_SINGLE(&آرگومانها, عدم, عدم, sa, sb, 0);
    }

#اگر_تعریف_شده SMP
  } وإلا {

    آرگومانها.n    = *N;
    info = GETRF_PARALLEL(&آرگومانها, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS > 0){
      آرگومانها.n    = *NRHS;
      GETRS_N_PARALLEL(&آرگومانها, عدم, عدم, sa, sb, 0);
    }
  }
#پایان_اگر

#اگر_تعریف_نشده PPC440
  blas_حافظه_free(میانگیر);
#پایان_اگر

  *Info = info;

  FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_انتها;

  أرجع;
}
