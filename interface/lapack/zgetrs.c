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
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XGETRS"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZGETRS"
#وإلا
#تعریف_کن ERROR_NAME "CGETRS"
#پایان_اگر

ساكن blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_SINGLE, GETRS_T_SINGLE, GETRS_R_SINGLE, GETRS_C_SINGLE,
};

#اگر_تعریف_شده SMP
ساكن blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  GETRS_N_PARALLEL, GETRS_T_PARALLEL, GETRS_R_PARALLEL, GETRS_C_PARALLEL,
};
#پایان_اگر

OPENBLAS_EXPORT
فراغ NAME(حرف *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
	    blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  حرف trans_arg = *TRANS;

  blas_arg_t آرگومانها;

  blasint info;
  صحيح trans;
  FLOAT *میانگیر;
#اگر_تعریف_شده PPC440
  خارجي
#پایان_اگر
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  آرگومانها.m    = *N;
  آرگومانها.n    = *NRHS;
  آرگومانها.a    = (فراغ *)a;
  آرگومانها.lda  = *ldA;
  آرگومانها.b    = (فراغ *)b;
  آرگومانها.ldb  = *ldB;
  آرگومانها.c    = (فراغ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 2;
  إذا (trans_arg == 'C') trans = 3;

  إذا (آرگومانها.ldb  < MAX(1, آرگومانها.m)) info = 8;
  إذا (آرگومانها.lda  < MAX(1, آرگومانها.m)) info = 5;
  إذا (آرگومانها.n    < 0) info = 3;
  إذا (آرگومانها.m    < 0) info = 2;
  إذا (trans     < 0) info = 1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    أرجع;
  }

  آرگومانها.alpha = عدم;
  آرگومانها.beta  = عدم;

  *Info = info;

  إذا (آرگومانها.m == 0 || آرگومانها.n == 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#اگر_تعریف_نشده PPC440
  میانگیر = (FLOAT *)blas_حافظه_alloc(1);

  sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);
#پایان_اگر

#اگر_تعریف_شده SMP
  إذا (آرگومانها.m < 64) {
    آرگومانها.nthreads = 1;
  } وإلا {
    آرگومانها.nthreads = رقم_cpu_avail(4);
  }

  إذا (آرگومانها.nthreads == 1) {
#پایان_اگر

    (getrs_single[trans])(&آرگومانها, عدم, عدم, sa, sb, 0);

#اگر_تعریف_شده SMP
  } وإلا {

    (getrs_parallel[trans])(&آرگومانها, عدم, عدم, sa, sb, 0);

  }
#پایان_اگر

#اگر_تعریف_نشده PPC440
  blas_حافظه_free(میانگیر);
#پایان_اگر

  FUNCTION_PROFILE_انتها(COMPSIZE * COMPSIZE, آرگومانها.m * آرگومانها.n, 2 * آرگومانها.m * آرگومانها.m * آرگومانها.n);

  IDEBUG_انتها;

  أرجع;

}
