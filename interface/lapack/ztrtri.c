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
#define ERROR_NAME "XTRTRI"
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZTRTRI"
#ଅନ୍ୟଥା
#define ERROR_NAME "CTRTRI"
#endif

ସ୍ଥିର blasint (*trtri_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UU_SINGLE, TRTRI_UN_SINGLE, TRTRI_LU_SINGLE, TRTRI_LN_SINGLE,
};

#ifdef SMP
ସ୍ଥିର blasint (*trtri_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UU_PARALLEL, TRTRI_UN_PARALLEL, TRTRI_LU_PARALLEL, TRTRI_LN_PARALLEL,
};
#endif

OPENBLAS_EXPORT
ଶୂନ୍ୟ NAME(ଅକ୍ଷର *UPLO, ଅକ୍ଷର *DIAG, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t args;

  blasint uplo_arg = *UPLO;
  blasint diag_arg = *DIAG;
  blasint uplo, diag;
  blasint info;
  FLOAT *buffer;
#ifdef PPC440
  ବାହ୍ୟ
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.n    = *N;
  args.a    = (ଶୂନ୍ୟ *)a;
  args.lda  = *ldA;

  TOUPPER(uplo_arg);
  TOUPPER(diag_arg);

  uplo = -1;
  ଯଦି (uplo_arg == 'U') uplo = 0;
  ଯଦି (uplo_arg == 'L') uplo = 1;
  diag = -1;
  ଯଦି (diag_arg == 'U') diag = 0;
  ଯଦି (diag_arg == 'N') diag = 1;

  info  = 0;
  ଯଦି (args.lda  < MAX(1,args.n)) info = 5;
  ଯଦି (args.n    < 0)             info = 3;
  ଯଦି (diag < 0)                  info = 2;
  ଯଦି (uplo < 0)                  info = 1;
  ଯଦି (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME) - 1);
    *Info = - info;
    ଫେରନ୍ତୁ;
  }

  *Info = 0;

  ଯଦି (args.n == 0) ଫେରନ୍ତୁ;

  ଯଦି (diag) {
    ଯଦି (AMIN_K(args.n, args.a, args.lda + 1) == ZERO) {
      *Info = IAMIN_K(args.n, args.a, args.lda + 1);
      ଫେରନ୍ତୁ;
    }
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  args.common = ଖାଲି;
  ଯଦି (args.n < 64) {
    args.nthreads = 1;
  } ଅନ୍ୟଥା {
    args.nthreads = ଅଙ୍କ_cpu_avail(4);
  }

  ଯଦି (args.nthreads == 1) {
#endif

    *Info = (trtri_single[(uplo << 1) | diag])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

#ifdef SMP
  } ଅନ୍ୟଥା {

    *Info = (trtri_parallel[(uplo << 1) | diag])(&args, ଖାଲି, ଖାଲି, sa, sb, 0);

  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(1, .5 * args.n * args.n,
		          2. * args.n * (1./3. + args.n * ( 1./2. + args.n * 1./6.))
		       +  6. * args.n * (1./3. + args.n * (-1./2. + args.n * 1./6.)));

  IDEBUG_END;

  ଫେରନ୍ତୁ;
}
