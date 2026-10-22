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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGETRS"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGETRS"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGETRS"
#ಆದರೆ_ಅಂತ್ಯ

ಸ್ಥಿರ blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_SINGLE, GETRS_T_SINGLE,
};

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
ಸ್ಥಿರ blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_PARALLEL, GETRS_T_PARALLEL,
};
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT
ಶೂನ್ಯ NAME(ಅಕ್ಷರ *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  ಅಕ್ಷರ trans_arg = *TRANS;

  blas_arg_t ವಾದಗಳು;

  blasint info;
  ಪೂರ್ಣಾಂಕ trans;
  FLOAT *ಬಫರ್;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ PPC440
  ಬಾಹ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  ವಾದಗಳು.m    = *N;
  ವಾದಗಳು.n    = *NRHS;
  ವಾದಗಳು.a    = (ಶೂನ್ಯ *)a;
  ವಾದಗಳು.lda  = *ldA;
  ವಾದಗಳು.b    = (ಶೂನ್ಯ *)b;
  ವಾದಗಳು.ldb  = *ldB;
  ವಾದಗಳು.c    = (ಶೂನ್ಯ *)ipiv;

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  ಆದರೆ (trans_arg == 'N') trans = 0;
  ಆದರೆ (trans_arg == 'T') trans = 1;
  ಆದರೆ (trans_arg == 'R') trans = 0;
  ಆದರೆ (trans_arg == 'C') trans = 1;

  ಆದರೆ (ವಾದಗಳು.ldb  < MAX(1, ವಾದಗಳು.m)) info = 8;
  ಆದರೆ (ವಾದಗಳು.lda  < MAX(1, ವಾದಗಳು.m)) info = 5;
  ಆದರೆ (ವಾದಗಳು.n    < 0) info = 3;
  ಆದರೆ (ವಾದಗಳು.m    < 0) info = 2;
  ಆದರೆ (trans     < 0) info = 1;

  ಆದರೆ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME) - 1);
    ಹಿಂತಿರುಗಿ;
  }

  ವಾದಗಳು.alpha = ಖಾಲಿ;
  ವಾದಗಳು.beta  = ಖಾಲಿ;

  *Info = info;

  ಆದರೆ (ವಾದಗಳು.m == 0 || ವಾದಗಳು.n == 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ PPC440
  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(1);

  sa = (FLOAT *)((BLASLONG)ಬಫರ್ + GEMM_ಸ್ಥಾನಾಂತರ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ವಾದಗಳು.common = ಖಾಲಿ;
  ಆದರೆ (ವಾದಗಳು.m < 64) {
    ವಾದಗಳು.nthreads = 1;
  } ಇಲ್ಲದಿದ್ದರೆ {
    ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(4);
  }

  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    (getrs_single[trans])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {
    (getrs_parallel[trans])(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);
  }
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ PPC440
  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, ವಾದಗಳು.m * ವಾದಗಳು.n, 2 * ವಾದಗಳು.m * ವಾದಗಳು.m * ವಾದಗಳು.n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;

}
