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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGESV"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGESV"
#وإلا
#define ERROR_NAME "SGESV"
#endif
#وإلا
#ifdef XDOUBLE
#define ERROR_NAME "XGESV"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZGESV"
#وإلا
#define ERROR_NAME "CGESV"
#endif
#endif

OPENBLAS_EXPORT
فراغ NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t حجوجي;

  blasint info;
  FLOAT *buffer;
#ifdef PPC440
  خارجي
#endif
  FLOAT *sa, *sb;

  بغ_DEBUG_NAME;

  حجوجي.m    = *N;
  حجوجي.n    = *NRHS;
  حجوجي.a    = (فراغ *)a;
  حجوجي.lda  = *ldA;
  حجوجي.b    = (فراغ *)b;
  حجوجي.ldb  = *ldB;
  حجوجي.c    = (فراغ *)ipiv;

  info  = 0;
  إذا (حجوجي.ldb < MAX(1,حجوجي.m)) info = 7;
  إذا (حجوجي.lda < MAX(1,حجوجي.m)) info = 4;
  إذا (حجوجي.n   < 0)             info = 2;
  إذا (حجوجي.m   < 0)             info = 1;

  إذا (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    أرجع;
  }

  حجوجي.alpha = عدم;
  حجوجي.beta  = عدم;

  *Info = 0;

  إذا (حجوجي.m == 0) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  حجوجي.common = عدم;

#إذا defined(_WIN64) && defined(_M_ARM64)
  #ifdef COMPLEX
    إذا (حجوجي.m * حجوجي.m <= 300) 
  #وإلا
    إذا (حجوجي.m * حجوجي.m <= 500) 
  #endif
      حجوجي.nthreads = 1;
  وإلا إذا (حجوجي.m * حجوجي.m <= 1000)
    حجوجي.nthreads = 4;
  وإلا
    حجوجي.nthreads = لمب_cpu_avail(4);
#وإلا
  #ifndef DOUBLE
    إذا (حجوجي.m * حجوجي.m < 40000)
  #وإلا
    إذا (حجوجي.m * حجوجي.m < 10000)
  #endif
      حجوجي.nthreads = 1;
    وإلا
      حجوجي.nthreads = لمب_cpu_avail(4);
#endif

  إذا (حجوجي.nthreads == 1) {
#endif

    حجوجي.n    = *N;
    info = GETRF_SINGLE(&حجوجي, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS >0){
      حجوجي.n    = *NRHS;
      GETRS_N_SINGLE(&حجوجي, عدم, عدم, sa, sb, 0);
    }

#ifdef SMP
  } وإلا {

    حجوجي.n    = *N;
    info = GETRF_PARALLEL(&حجوجي, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS > 0){
      حجوجي.n    = *NRHS;
      GETRS_N_PARALLEL(&حجوجي, عدم, عدم, sa, sb, 0);
    }
  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  *Info = info;

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_END;

  أرجع;
}
