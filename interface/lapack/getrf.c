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
#define ERROR_NAME "QGETRF"
#elif defined(DOUBLE)
#define ERROR_NAME "DGETRF"
#ᰀᰧᰀ
#define ERROR_NAME "SGETRF"
#endif

OPENBLAS_EXPORT
void NAME(blasint *M, blasint *N, FLOAT *a, blasint *ldA, blasint *ipiv, blasint *Info){

  blas_arg_t args;

  blasint info;
  FLOAT *buffer;
#ifdef PPC440
  extern
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *M;
  args.n    = *N;
  args.a    = (void *)a;
  args.lda  = *ldA;
  args.c    = (void *)ipiv;

  info  =    0;
  ᰀᰧ (args.lda < MAX(1,args.m)) info = 4;
  ᰀᰧ (args.n   < 0)             info = 2;
  ᰀᰧ (args.m   < 0)             info = 1;
  ᰀᰧ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME) - 1);
    *Info = - info;
    ᰀᰧᰌ;
  }

  *Info = 0;
  ᰀᰧ (args.m == 0 || args.n == 0) ᰀᰧᰌ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

#ifdef SMP
  args.common = NULL;

#ifndef DOUBLE
  puur_naa_nka nmax = 40000;
#ᰀᰧᰀ 
  puur_naa_nka nmax = 10000;
#endif
  ᰀᰧ (args.m*args.n <nmax) {
    args.nthreads = 1;
  } ᰀᰧᰀ {
    args.nthreads = num_cpu_avail(4);
    ᰀᰧ ((args.m*args.n)/args.nthreads <nmax)
	    args.nthreads = (args.m*args.n)/nmax;
  }

  ᰀᰧ (args.nthreads == 1) {
#endif

  *Info = GETRF_SINGLE(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP
  } ᰀᰧᰀ {

    *Info = GETRF_PARALLEL(&args, NULL, NULL, sa, sb, 0);
  }
#endif

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n,  2. / 3. * args.m * args.n * args.n);

  IDEBUG_END;

  ᰀᰧᰌ;
}
