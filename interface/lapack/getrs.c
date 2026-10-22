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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGETRS"
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGETRS"
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGETRS"
#endif

സ്ഥിതം blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_SINGLE, GETRS_T_SINGLE,
};

#ifdef SMP
സ്ഥിതം blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_PARALLEL, GETRS_T_PARALLEL,
};
#endif

OPENBLAS_EXPORT
ശൂന്യം NAME(അക്ഷരം *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  അക്ഷരം trans_arg = *TRANS;

  blas_arg_t args;

  blasint info;
  പൂർണ്ണസംഖ്യ trans;
  FLOAT *buffer;
#ifdef PPC440
  ബാഹ്യ
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *N;
  args.n    = *NRHS;
  args.a    = (ശൂന്യം *)a;
  args.lda  = *ldA;
  args.b    = (ശൂന്യം *)b;
  args.ldb  = *ldB;
  args.c    = (ശൂന്യം *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  എങ്കിൽ (trans_arg == 'N') trans = 0;
  എങ്കിൽ (trans_arg == 'T') trans = 1;
  എങ്കിൽ (trans_arg == 'R') trans = 0;
  എങ്കിൽ (trans_arg == 'C') trans = 1;

  എങ്കിൽ (args.ldb  < MAX(1, args.m)) info = 8;
  എങ്കിൽ (args.lda  < MAX(1, args.m)) info = 5;
  എങ്കിൽ (args.n    < 0) info = 3;
  എങ്കിൽ (args.m    < 0) info = 2;
  എങ്കിൽ (trans     < 0) info = 1;

  എങ്കിൽ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME) - 1);
    തിരികെ;
  }

  args.alpha = ശൂന്യ;
  args.beta  = ശൂന്യ;

  *Info = info;

  എങ്കിൽ (args.m == 0 || args.n == 0) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  args.common = ശൂന്യ;
  എങ്കിൽ (args.m < 64) {
    args.nthreads = 1;
  } അല്ലെങ്കിൽ {
    args.nthreads = അക്കം_cpu_avail(4);
  }

  എങ്കിൽ (args.nthreads == 1) {
#endif

    (getrs_single[trans])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);

#ifdef SMP
  } അല്ലെങ്കിൽ {
    (getrs_parallel[trans])(&args, ശൂന്യ, ശൂന്യ, sa, sb, 0);
  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n, 2 * args.m * args.m * args.n);

  IDEBUG_END;

  തിരികെ;

}
