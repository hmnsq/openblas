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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XGETRS"
#elif defined(DOUBLE)
#define ERROR_NAME "ZGETRS"
#else
#define ERROR_NAME "CGETRS"
#endif

static blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#ifdef SMP
static blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#endif

OPENBLAS_EXPORT
void NAME(char *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  char trans_arg = *TRANS;

  blas_arg_t मगा_कातविरान;

  blasint info;
  int trans;
  FLOAT *buffer;
#ifdef PPC440
  extern
#endif
  FLOAT *sa, *sb;

  इलिम्बाग_DEBUG_NAME;

  मगा_कातविरान.m    = *N;
  मगा_कातविरान.n    = *NRHS;
  मगा_कातविरान.a    = (void *)a;
  मगा_कातविरान.lda  = *ldA;
  मगा_कातविरान.b    = (void *)b;
  मगा_कातविरान.ldb  = *ldB;
  मगा_कातविरान.c    = (void *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  if (trans_arg == 'N') trans = 0;
  if (trans_arg == 'T') trans = 1;
  if (trans_arg == 'R') trans = 2;
  if (trans_arg == 'C') trans = 3;

  if (मगा_कातविरान.ldb  < MAX(1, मगा_कातविरान.m)) info = 8;
  if (मगा_कातविरान.lda  < MAX(1, मगा_कातविरान.m)) info = 5;
  if (मगा_कातविरान.n    < 0) info = 3;
  if (मगा_कातविरान.m    < 0) info = 2;
  if (trans     < 0) info = 1;

  if (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME) - 1);
    return;
  }

  मगा_कातविरान.alpha = NULL;
  मगा_कातविरान.beta  = NULL;

  *Info = info;

  if (मगा_कातविरान.m == 0 || मगा_कातविरान.n == 0) return;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  if (मगा_कातविरान.m < 64) {
    मगा_कातविरान.nthreads = 1;
  } else {
    मगा_कातविरान.nthreads = नुमेरो_cpu_avail(4);
  }

  if (मगा_कातविरान.nthreads == 1) {
#endif

    (getrs_single[trans])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } else {

    (getrs_parallel[trans])(&मगा_कातविरान, NULL, NULL, sa, sb, 0);

  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, मगा_कातविरान.m * मगा_कातविरान.n, 2 * मगा_कातविरान.m * मगा_कातविरान.m * मगा_कातविरान.n);

  IDEBUG_END;

  return;

}
