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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGESV"
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGESV"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGESV"
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGESV"
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGESV"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGESV"
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT
ਖਾਲੀ NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t ਦਲੀਲਾਂ;

  blasint info;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ PPC440
  ਬਾਹਰੀ
#ਜੇ_ਅੰਤ
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  ਦਲੀਲਾਂ.m    = *N;
  ਦਲੀਲਾਂ.n    = *NRHS;
  ਦਲੀਲਾਂ.a    = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.lda  = *ldA;
  ਦਲੀਲਾਂ.b    = (ਖਾਲੀ *)b;
  ਦਲੀਲਾਂ.ldb  = *ldB;
  ਦਲੀਲਾਂ.c    = (ਖਾਲੀ *)ipiv;

  info  = 0;
  ਜੇ (ਦਲੀਲਾਂ.ldb < MAX(1,ਦਲੀਲਾਂ.m)) info = 7;
  ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,ਦਲੀਲਾਂ.m)) info = 4;
  ਜੇ (ਦਲੀਲਾਂ.n   < 0)             info = 2;
  ਜੇ (ਦਲੀਲਾਂ.m   < 0)             info = 1;

  ਜੇ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME) - 1);
    *Info = - info;
    ਵਾਪਸ;
  }

  ਦਲੀਲਾਂ.alpha = ਨਿਰਾ;
  ਦਲੀਲਾਂ.beta  = ਨਿਰਾ;

  *Info = 0;

  ਜੇ (ਦਲੀਲਾਂ.m == 0) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ PPC440
  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)ਬਫਰ + GEMM_ਵਿੱਥ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਦਲੀਲਾਂ.common = ਨਿਰਾ;

#ਜੇ defined(_WIN64) && defined(_M_ARM64)
  #ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m <= 300) 
  #ਨਹੀਂਤਾਂ
    ਜੇ (ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m <= 500) 
  #ਜੇ_ਅੰਤ
      ਦਲੀਲਾਂ.nthreads = 1;
  ਨਹੀਂਤਾਂ ਜੇ (ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m <= 1000)
    ਦਲੀਲਾਂ.nthreads = 4;
  ਨਹੀਂਤਾਂ
    ਦਲੀਲਾਂ.nthreads = ਅੰਕ_cpu_avail(4);
#ਨਹੀਂਤਾਂ
  #ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ DOUBLE
    ਜੇ (ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m < 40000)
  #ਨਹੀਂਤਾਂ
    ਜੇ (ਦਲੀਲਾਂ.m * ਦਲੀਲਾਂ.m < 10000)
  #ਜੇ_ਅੰਤ
      ਦਲੀਲਾਂ.nthreads = 1;
    ਨਹੀਂਤਾਂ
      ਦਲੀਲਾਂ.nthreads = ਅੰਕ_cpu_avail(4);
#ਜੇ_ਅੰਤ

  ਜੇ (ਦਲੀਲਾਂ.nthreads == 1) {
#ਜੇ_ਅੰਤ

    ਦਲੀਲਾਂ.n    = *N;
    info = GETRF_SINGLE(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

    ਜੇ (info == 0 && *NRHS >0){
      ਦਲੀਲਾਂ.n    = *NRHS;
      GETRS_N_SINGLE(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);
    }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

    ਦਲੀਲਾਂ.n    = *N;
    info = GETRF_PARALLEL(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

    ਜੇ (info == 0 && *NRHS > 0){
      ਦਲੀਲਾਂ.n    = *NRHS;
      GETRS_N_PARALLEL(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);
    }
  }
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ PPC440
  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);
#ਜੇ_ਅੰਤ

  *Info = info;

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
