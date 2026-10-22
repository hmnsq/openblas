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

#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGETRS"
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGETRS"
#وإلا
#وصف_ڪريو ERROR_NAME "CGETRS"
#جيڪڏهن_پورو

ساكن blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#جيڪڏهن_وصف_ٿيل SMP
ساكن blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#جيڪڏهن_پورو

OPENBLAS_EXPORT
فراغ NAME(حرف *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  حرف trans_arg = *TRANS;

  blas_arg_t دليلون;

  blasint info;
  صحيح trans;
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

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (دليلون.ldb  < MAX(1, دليلون.m)) info = 8;
  إذا (دليلون.lda  < MAX(1, دليلون.m)) info = 5;
  إذا (دليلون.n    < 0) info = 3;
  إذا (دليلون.m    < 0) info = 2;
  إذا (trans     < 0) info = 1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    أرجع;
  }

  دليلون.alpha = عدم;
  دليلون.beta  = عدم;

  *Info = info;

  إذا (دليلون.m == 0 || دليلون.n == 0) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

#جيڪڏهن_وصف_نه_ٿيل PPC440
  بفر = (FLOAT *)blas_يادگيري_alloc(1);

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
  إذا (دليلون.m < 64) {
    دليلون.nthreads = 1;
  } وإلا {
    دليلون.nthreads = انگ_cpu_avail(4);
  }

  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    (getrs_single[trans])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

    (getrs_parallel[trans])(&دليلون, عدم, عدم, sa, sb, 0);

  }
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل PPC440
  blas_يادگيري_free(بفر);
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE, دليلون.m * دليلون.n, 2 * دليلون.m * دليلون.m * دليلون.n);

  IDEBUG_پڇاڙي;

  أرجع;

}
