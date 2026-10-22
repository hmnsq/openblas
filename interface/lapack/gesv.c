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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QGESV"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DGESV"
#وإلا
#وصف_ڪريو ERROR_NAME "SGESV"
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGESV"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGESV"
#وإلا
#وصف_ڪريو ERROR_NAME "CGESV"
#جيڪڏهن_پورو
#جيڪڏهن_پورو

OPENBLAS_EXPORT
فراغ NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t دليلون;

  blasint info;
  FLOAT *بفر;
#جيڪڏهن_وصف_ٿيل PPC440
  خارجي
#جيڪڏهن_پورو
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  دليلون.m    = *N;
  دليلون.n    = *NRHS;
  دليلون.a    = (فراغ *)a;
  دليلون.lda  = *ldA;
  دليلون.b    = (فراغ *)b;
  دليلون.ldb  = *ldB;
  دليلون.c    = (فراغ *)ipiv;

  info  = 0;
  إذا (دليلون.ldb < MAX(1,دليلون.m)) info = 7;
  إذا (دليلون.lda < MAX(1,دليلون.m)) info = 4;
  إذا (دليلون.n   < 0)             info = 2;
  إذا (دليلون.m   < 0)             info = 1;

  إذا (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    أرجع;
  }

  دليلون.alpha = عدم;
  دليلون.beta  = عدم;

  *Info = 0;

  إذا (دليلون.m == 0) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

#جيڪڏهن_وصف_نه_ٿيل PPC440
  بفر = (FLOAT *)blas_يادگيري_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
  دليلون.common = عدم;

#إذا defined(_WIN64) && defined(_M_ARM64)
  #جيڪڏهن_وصف_ٿيل COMPLEX
    إذا (دليلون.m * دليلون.m <= 300) 
  #وإلا
    إذا (دليلون.m * دليلون.m <= 500) 
  #جيڪڏهن_پورو
      دليلون.nthreads = 1;
  وإلا إذا (دليلون.m * دليلون.m <= 1000)
    دليلون.nthreads = 4;
  وإلا
    دليلون.nthreads = انگ_cpu_avail(4);
#وإلا
  #جيڪڏهن_وصف_نه_ٿيل DOUBLE
    إذا (دليلون.m * دليلون.m < 40000)
  #وإلا
    إذا (دليلون.m * دليلون.m < 10000)
  #جيڪڏهن_پورو
      دليلون.nthreads = 1;
    وإلا
      دليلون.nthreads = انگ_cpu_avail(4);
#جيڪڏهن_پورو

  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    دليلون.n    = *N;
    info = GETRF_SINGLE(&دليلون, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS >0){
      دليلون.n    = *NRHS;
      GETRS_N_SINGLE(&دليلون, عدم, عدم, sa, sb, 0);
    }

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

    دليلون.n    = *N;
    info = GETRF_PARALLEL(&دليلون, عدم, عدم, sa, sb, 0);

    إذا (info == 0 && *NRHS > 0){
      دليلون.n    = *NRHS;
      GETRS_N_PARALLEL(&دليلون, عدم, عدم, sa, sb, 0);
    }
  }
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل PPC440
  blas_يادگيري_free(بفر);
#جيڪڏهن_پورو

  *Info = info;

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_پڇاڙي;

  أرجع;
}
