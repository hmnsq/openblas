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
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGESV"
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGESV"
#ورنہ
#تعریف_کرو ERROR_NAME "SGESV"
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGESV"
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGESV"
#ورنہ
#تعریف_کرو ERROR_NAME "CGESV"
#اگر_ختم
#اگر_ختم

OPENBLAS_EXPORT
خالی NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t دلائل;

  blasint info;
  FLOAT *بفر;
#اگر_تعریف_شدہ PPC440
  خارجی
#اگر_ختم
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  دلائل.m    = *N;
  دلائل.n    = *NRHS;
  دلائل.a    = (خالی *)a;
  دلائل.lda  = *ldA;
  دلائل.b    = (خالی *)b;
  دلائل.ldb  = *ldB;
  دلائل.c    = (خالی *)ipiv;

  info  = 0;
  اگر (دلائل.ldb < MAX(1,دلائل.m)) info = 7;
  اگر (دلائل.lda < MAX(1,دلائل.m)) info = 4;
  اگر (دلائل.n   < 0)             info = 2;
  اگر (دلائل.m   < 0)             info = 1;

  اگر (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    واپس;
  }

  دلائل.alpha = تہی;
  دلائل.beta  = تہی;

  *Info = 0;

  اگر (دلائل.m == 0) واپس;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

#اگر_غیر_تعریف_شدہ PPC440
  بفر = (FLOAT *)blas_حافظہ_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);
#اگر_ختم

#اگر_تعریف_شدہ SMP
  دلائل.common = تہی;

#اگر defined(_WIN64) && defined(_M_ARM64)
  #اگر_تعریف_شدہ COMPLEX
    اگر (دلائل.m * دلائل.m <= 300) 
  #ورنہ
    اگر (دلائل.m * دلائل.m <= 500) 
  #اگر_ختم
      دلائل.nthreads = 1;
  ورنہ اگر (دلائل.m * دلائل.m <= 1000)
    دلائل.nthreads = 4;
  ورنہ
    دلائل.nthreads = عدد_cpu_avail(4);
#ورنہ
  #اگر_غیر_تعریف_شدہ DOUBLE
    اگر (دلائل.m * دلائل.m < 40000)
  #ورنہ
    اگر (دلائل.m * دلائل.m < 10000)
  #اگر_ختم
      دلائل.nthreads = 1;
    ورنہ
      دلائل.nthreads = عدد_cpu_avail(4);
#اگر_ختم

  اگر (دلائل.nthreads == 1) {
#اگر_ختم

    دلائل.n    = *N;
    info = GETRF_SINGLE(&دلائل, تہی, تہی, sa, sb, 0);

    اگر (info == 0 && *NRHS >0){
      دلائل.n    = *NRHS;
      GETRS_N_SINGLE(&دلائل, تہی, تہی, sa, sb, 0);
    }

#اگر_تعریف_شدہ SMP
  } ورنہ {

    دلائل.n    = *N;
    info = GETRF_PARALLEL(&دلائل, تہی, تہی, sa, sb, 0);

    اگر (info == 0 && *NRHS > 0){
      دلائل.n    = *NRHS;
      GETRS_N_PARALLEL(&دلائل, تہی, تہی, sa, sb, 0);
    }
  }
#اگر_ختم

#اگر_غیر_تعریف_شدہ PPC440
  blas_حافظہ_free(بفر);
#اگر_ختم

  *Info = info;

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_انتہا;

  واپس;
}
