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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGETRS"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGETRS"
#ଅନ୍ୟଥା
#define ERROR_NAME "CGETRS"
#endif

ସ୍ଥିର blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#ifdef SMP
ସ୍ଥିର blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#endif

OPENBLAS_EXPORT
ଶୂନ୍ୟ NAME(ଅକ୍ଷର *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  ଅକ୍ଷର trans_arg = *TRANS;

  blas_arg_t args;

  blasint info;
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା trans;
  FLOAT *buffer;
#ifdef PPC440
  ବାହ୍ୟ
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *N;
  args.n    = *NRHS;
  args.a    = (ଶୂନ୍ୟ *)a;
  args.lda  = *ldA;
  args.b    = (ଶୂନ୍ୟ *)b;
  args.ldb  = *ldB;
  args.c    = (ଶୂନ୍ୟ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  ଯଦି (trans_arg == 'N') trans = 0;
  ଯଦି (trans_arg == 'T') trans = 1;
  ଯଦି (trans_arg == 'R') trans = 2;
  ଯଦି (trans_arg == 'C') trans = 3;

  ଯଦି (args.ldb  < MAX(1, args.m)) info = 8;
  ଯଦି (args.lda  < MAX(1, args.m)) info = 5;
  ଯଦି (args.n    < 0) info = 3;
  ଯଦି (args.m    < 0) info = 2;
  ଯଦି (trans     < 0) info = 1;

  ଯଦି (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME) - 1);
    ଫେରନ୍ତୁ;
  }

  args.alpha = ଖାଲି;
  args.beta  = ଖାଲି;

  *Info = info;

  ଯଦି (args.m == 0 || args.n == 0) ଫେରନ୍ତୁ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  ଯଦି (args.m < 64) {
    args.nthreads = 1;
  } ଅନ୍ୟଥା {
    args.nthreads = ଅଙ୍କ_cpu_avail(4);
  }

  ଯଦି (args.nthreads == 1) {
#endif

    (getrs_single[trans])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP
  } ଅନ୍ୟଥା {

    (getrs_parallel[trans])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n, 2 * args.m * args.m * args.n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;

}
