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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGETRS"
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGETRS"
#وإلا
#تعریف_کرو ERROR_NAME "CGETRS"
#جے_انت

ساكن blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#جے_تعریف_شدہ SMP
ساكن blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#جے_انت

OPENBLAS_EXPORT
فراغ NAME(حرف *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  حرف trans_arg = *TRANS;

  blas_arg_t دلیلاں;

  blasint info;
  صحيح trans;
  FLOAT *بفر;
#جے_تعریف_شدہ PPC440
  خارجي
#جے_انت
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  دلیلاں.m    = *N;
  دلیلاں.n    = *NRHS;
  دلیلاں.a    = (فراغ *)a;
  دلیلاں.lda  = *ldA;
  دلیلاں.b    = (فراغ *)b;
  دلیلاں.ldb  = *ldB;
  دلیلاں.c    = (فراغ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (دلیلاں.ldb  < MAX(1, دلیلاں.m)) info = 8;
  إذا (دلیلاں.lda  < MAX(1, دلیلاں.m)) info = 5;
  إذا (دلیلاں.n    < 0) info = 3;
  إذا (دلیلاں.m    < 0) info = 2;
  إذا (trans     < 0) info = 1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    أرجع;
  }

  دلیلاں.alpha = عدم;
  دلیلاں.beta  = عدم;

  *Info = info;

  إذا (دلیلاں.m == 0 || دلیلاں.n == 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#جے_بنا_تعریف PPC440
  بفر = (FLOAT *)blas_حافظہ_alloc(1);

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وتھ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ناپ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وتھ_B);
#جے_انت

#جے_تعریف_شدہ SMP
  إذا (دلیلاں.m < 64) {
    دلیلاں.nthreads = 1;
  } وإلا {
    دلیلاں.nthreads = عدد_cpu_avail(4);
  }

  إذا (دلیلاں.nthreads == 1) {
#جے_انت

    (getrs_single[trans])(&دلیلاں, عدم, عدم, sa, sb, 0);

#جے_تعریف_شدہ SMP
  } وإلا {

    (getrs_parallel[trans])(&دلیلاں, عدم, عدم, sa, sb, 0);

  }
#جے_انت

#جے_بنا_تعریف PPC440
  blas_حافظہ_free(بفر);
#جے_انت

  FUNCTION_PROFILE_سرا(COMPSIZE * COMPSIZE, دلیلاں.m * دلیلاں.n, 2 * دلیلاں.m * دلیلاں.m * دلیلاں.n);

  IDEBUG_سرا;

  أرجع;

}
