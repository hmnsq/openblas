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

#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTRTRS"
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTRTRS"
#ورنہ
#تعریف_کرو ERROR_NAME "STRTRS"
#اگر_ختم

ساکن blasint (*trtrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
    TRTRS_UNU_SINGLE, TRTRS_UNN_SINGLE, TRTRS_UTU_SINGLE, TRTRS_UTN_SINGLE, TRTRS_LNU_SINGLE, TRTRS_LNN_SINGLE, TRTRS_LTU_SINGLE, TRTRS_LTN_SINGLE,
};

#اگر_تعریف_شدہ SMP
ساکن blasint (*trtrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
    TRTRS_UNU_PARALLEL, TRTRS_UNN_PARALLEL, TRTRS_UTU_PARALLEL, TRTRS_UTN_PARALLEL, TRTRS_LNU_PARALLEL, TRTRS_LNN_PARALLEL, TRTRS_LTU_PARALLEL, TRTRS_LTN_PARALLEL,
};
#اگر_ختم

OPENBLAS_EXPORT
خالی NAME(حرف *UPLO, حرف* TRANS, حرف* DIAG, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  FLOAT *b, blasint *ldB, blasint *Info){

    حرف uplo_arg = *UPLO;
    حرف trans_arg = *TRANS;
    حرف diag_arg = *DIAG;

  blas_arg_t دلائل;

  blasint info;
  صحیح uplo, trans, diag;
  FLOAT *بفر;
#اگر_تعریف_شدہ PPC440
  خارجی
#اگر_ختم
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  دلائل.m    = *N;
  دلائل.n    = *NRHS;
  دلائل.a    = (خالی *)a;
  دلائل.lda  = *ldA;
  دلائل.b    = (خالی *)b;
  دلائل.ldb  = *ldB;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'R') trans = 0;
  اگر (trans_arg == 'C') trans = 1;

  TOUPPER(uplo_arg);
  uplo = -1;
  اگر (uplo_arg == 'U') uplo = 0;
  اگر (uplo_arg == 'L') uplo = 1;

  TOUPPER(diag_arg);
  diag = -1;
  اگر (diag_arg == 'U') diag = 0;
  اگر (diag_arg == 'N') diag = 1;

  اگر (دلائل.ldb  < MAX(1, دلائل.m)) info = 9;
  اگر (دلائل.lda  < MAX(1, دلائل.m)) info = 7;
  اگر (دلائل.n    < 0) info = 5;
  اگر (دلائل.m    < 0) info = 4;
  اگر (trans     < 0) info = 2;
  اگر (uplo      < 0) info = 1;
  اگر (diag      < 0) info = 3;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    واپس;
  }

  دلائل.alpha = تہی;
  دلائل.beta  = تہی;

  *Info = 0;

  اگر (دلائل.m == 0) واپس;

  اگر (diag) {
    اگر (AMIN_K(دلائل.m, دلائل.a, دلائل.lda + 1) == ZERO) {
      *Info = IAMIN_K(دلائل.m, دلائل.a, دلائل.lda + 1);
      واپس;
    }
  }


  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

#اگر_غیر_تعریف_شدہ PPC440
  بفر = (FLOAT *)blas_حافظہ_alloc(1);

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);
#اگر_ختم

#اگر_تعریف_شدہ SMP
  دلائل.common = تہی;
  اگر (دلائل.m < 64) {
    دلائل.nthreads = 1;
  } ورنہ {
    دلائل.nthreads = عدد_cpu_avail(4);
  }

  اگر (دلائل.nthreads == 1) {
#اگر_ختم

      (trtrs_single[(uplo << 2) | (trans << 1) | diag])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP
  } ورنہ {
    (trtrs_parallel[(uplo << 2) | (trans << 1) | diag])(&دلائل, تہی, تہی, sa, sb, 0);
  }
#اگر_ختم

#اگر_غیر_تعریف_شدہ PPC440
  blas_حافظہ_free(بفر);
#اگر_ختم

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE, دلائل.m * دلائل.n, 2 * دلائل.m * دلائل.m * دلائل.n);

  IDEBUG_انتہا;

  واپس;

}
