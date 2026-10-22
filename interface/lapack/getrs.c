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

#ifdef XDOUBLE
#define ERROR_NAME "QGETRS"
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGETRS"
#ܐܚܪܢܐ
#define ERROR_NAME "SGETRS"
#endif

ساكن blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_SINGLE, GETRS_T_SINGLE,
};

#ifdef SMP
ساكن blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_PARALLEL, GETRS_T_PARALLEL,
};
#endif

OPENBLAS_EXPORT
فراغ NAME(حرف *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  حرف trans_arg = *TRANS;

  blas_arg_t args;

  blasint info;
  صحيح trans;
  FLOAT *buffer;
#ifdef PPC440
  خارجي
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *N;
  args.n    = *NRHS;
  args.a    = (فراغ *)a;
  args.lda  = *ldA;
  args.b    = (فراغ *)b;
  args.ldb  = *ldB;
  args.c    = (فراغ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  ܐܢ (trans_arg == 'N') trans = 0;
  ܐܢ (trans_arg == 'T') trans = 1;
  ܐܢ (trans_arg == 'R') trans = 0;
  ܐܢ (trans_arg == 'C') trans = 1;

  ܐܢ (args.ldb  < MAX(1, args.m)) info = 8;
  ܐܢ (args.lda  < MAX(1, args.m)) info = 5;
  ܐܢ (args.n    < 0) info = 3;
  ܐܢ (args.m    < 0) info = 2;
  ܐܢ (trans     < 0) info = 1;

  ܐܢ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    ܦܢܐ;
  }

  args.alpha = عدم;
  args.beta  = عدم;

  *Info = info;

  ܐܢ (args.m == 0 || args.n == 0) ܦܢܐ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  args.common = عدم;
  ܐܢ (args.m < 64) {
    args.nthreads = 1;
  } ܐܚܪܢܐ {
    args.nthreads = num_cpu_avail(4);
  }

  ܐܢ (args.nthreads == 1) {
#endif

    (getrs_single[trans])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } ܐܚܪܢܐ {
    (getrs_parallel[trans])(&args, عدم, عدم, sa, sb, 0);
  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n, 2 * args.m * args.m * args.n);

  IDEBUG_END;

  ܦܢܐ;

}
