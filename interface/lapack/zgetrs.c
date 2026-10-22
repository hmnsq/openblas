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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGETRS"
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGETRS"
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGETRS"
#endif

ಸ್ಥಿರ blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#ifdef SMP
ಸ್ಥಿರ blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#endif

OPENBLAS_EXPORT
ಶೂನ್ಯ NAME(ವರ್ಣ *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  ವರ್ಣ trans_arg = *TRANS;

  blas_arg_t args;

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ trans;
  FLOAT *buffer;
#ifdef PPC440
  ಬಾಹ್ಯ
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *N;
  args.n    = *NRHS;
  args.a    = (ಶೂನ್ಯ *)a;
  args.lda  = *ldA;
  args.b    = (ಶೂನ್ಯ *)b;
  args.ldb  = *ldB;
  args.c    = (ಶೂನ್ಯ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'R') trans = 2;
  ಯದಿ (trans_arg == 'C') trans = 3;

  ಯದಿ (args.ldb  < MAX(1, args.m)) info = 8;
  ಯದಿ (args.lda  < MAX(1, args.m)) info = 5;
  ಯದಿ (args.n    < 0) info = 3;
  ಯದಿ (args.m    < 0) info = 2;
  ಯದಿ (trans     < 0) info = 1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME) - 1);
    ಪ್ರತ್ಯಾಗಮ;
  }

  args.alpha = ರಿಕ್ತ;
  args.beta  = ರಿಕ್ತ;

  *Info = info;

  ಯದಿ (args.m == 0 || args.n == 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  ಯದಿ (args.m < 64) {
    args.nthreads = 1;
  } ಅನ್ಯಥಾ {
    args.nthreads = ಅಂಕೆ_cpu_avail(4);
  }

  ಯದಿ (args.nthreads == 1) {
#endif

    (getrs_single[trans])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP
  } ಅನ್ಯಥಾ {

    (getrs_parallel[trans])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n, 2 * args.m * args.m * args.n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;

}
