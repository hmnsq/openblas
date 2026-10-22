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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QPOTRI"
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DPOTRI"
#નહીંતર
#define ERROR_NAME "SPOTRI"
#endif

સ્થિર blasint (*trtri_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UN_SINGLE, TRTRI_LN_SINGLE,
};

સ્થિર blasint (*lauum_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  LAUUM_U_SINGLE, LAUUM_L_SINGLE,
};

#ifdef SMP
સ્થિર blasint (*trtri_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UN_PARALLEL, TRTRI_LN_PARALLEL,
};

સ્થિર blasint (*lauum_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  LAUUM_U_PARALLEL, LAUUM_L_PARALLEL,
};
#endif

OPENBLAS_EXPORT
શૂન્ય NAME(અક્ષર *UPLO, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t args;

  blasint uplo_arg = *UPLO;
  blasint uplo;
  blasint info;
  FLOAT *buffer;
#ifdef PPC440
  બાહ્ય
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.n    = *N;
  args.a    = (શૂન્ય *)a;
  args.lda  = *ldA;

  TOUPPER(uplo_arg);

  uplo = -1;
  જો (uplo_arg == 'U') uplo = 0;
  જો (uplo_arg == 'L') uplo = 1;

  info  = 0;
  જો (args.lda  < MAX(1,args.n)) info = 4;
  જો (args.n    < 0)             info = 2;
  જો (uplo < 0)                  info = 1;

  જો (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME) - 1);
    *Info = - info;
    પાછા;
  }

  *Info = 0;

  જો (args.n == 0) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  args.common = ખાલી;
  જો (args.n < 180)
    args.nthreads = 1;
  નહીંતર
  args.nthreads = અંક_cpu_avail(4);

  જો (args.nthreads == 1) {
#endif

  info = (trtri_single[uplo])(&args, ખાલી, ખાલી, sa, sb, 0);

  જો (!info) {
    info = (lauum_single[uplo])(&args, ખાલી, ખાલી, sa, sb, 0);
  }

  *Info = info;

#ifdef SMP
  } નહીંતર {
    info = (trtri_parallel[uplo])(&args, ખાલી, ખાલી, sa, sb, 0);

    જો (!info) {
      info = (lauum_parallel[uplo])(&args, ખાલી, ખાલી, sa, sb, 0);
    }

    *Info = info;
  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, .5 * args.n * args.n,
		          args.n * (1./3. + args.n * ( 1./2. + args.n * 1./6.))
		       +  args.n * (1./3. + args.n * (-1./2. + args.n * 1./6.)));

  IDEBUG_END;

  પાછા;
}
