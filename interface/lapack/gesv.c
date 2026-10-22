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

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGESV"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGESV"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGESV"
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGESV"
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGESV"
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGESV"
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT
ಶೂನ್ಯ NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t ವಾದಗಳು;

  blasint info;
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

  info  = 0;
  ಆದರೆ (ವಾದಗಳು.ldb < MAX(1,ವಾದಗಳು.m)) info = 7;
  ಆದರೆ (ವಾದಗಳು.lda < MAX(1,ವಾದಗಳು.m)) info = 4;
  ಆದರೆ (ವಾದಗಳು.n   < 0)             info = 2;
  ಆದರೆ (ವಾದಗಳು.m   < 0)             info = 1;

  ಆದರೆ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME) - 1);
    *Info = - info;
    ಹಿಂತಿರುಗಿ;
  }

  ವಾದಗಳು.alpha = ಖಾಲಿ;
  ವಾದಗಳು.beta  = ಖಾಲಿ;

  *Info = 0;

  ಆದರೆ (ವಾದಗಳು.m == 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ PPC440
  ಬಫರ್ = (FLOAT *)blas_ಸ್ಮರಣೆ_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)ಬಫರ್ + GEMM_ಸ್ಥಾನಾಂತರ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ಗಾತ್ರ_2 + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ಸ್ಥಾನಾಂತರ_B);
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ವಾದಗಳು.common = ಖಾಲಿ;

#ಆದರೆ defined(_WIN64) && defined(_M_ARM64)
  #ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ COMPLEX
    ಆದರೆ (ವಾದಗಳು.m * ವಾದಗಳು.m <= 300) 
  #ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (ವಾದಗಳು.m * ವಾದಗಳು.m <= 500) 
  #ಆದರೆ_ಅಂತ್ಯ
      ವಾದಗಳು.nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (ವಾದಗಳು.m * ವಾದಗಳು.m <= 1000)
    ವಾದಗಳು.nthreads = 4;
  ಇಲ್ಲದಿದ್ದರೆ
    ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(4);
#ಇಲ್ಲದಿದ್ದರೆ
  #ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ DOUBLE
    ಆದರೆ (ವಾದಗಳು.m * ವಾದಗಳು.m < 40000)
  #ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (ವಾದಗಳು.m * ವಾದಗಳು.m < 10000)
  #ಆದರೆ_ಅಂತ್ಯ
      ವಾದಗಳು.nthreads = 1;
    ಇಲ್ಲದಿದ್ದರೆ
      ವಾದಗಳು.nthreads = ಅಂಕೆ_cpu_avail(4);
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (ವಾದಗಳು.nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

    ವಾದಗಳು.n    = *N;
    info = GETRF_SINGLE(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

    ಆದರೆ (info == 0 && *NRHS >0){
      ವಾದಗಳು.n    = *NRHS;
      GETRS_N_SINGLE(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);
    }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

    ವಾದಗಳು.n    = *N;
    info = GETRF_PARALLEL(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);

    ಆದರೆ (info == 0 && *NRHS > 0){
      ವಾದಗಳು.n    = *NRHS;
      GETRS_N_PARALLEL(&ವಾದಗಳು, ಖಾಲಿ, ಖಾಲಿ, sa, sb, 0);
    }
  }
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ PPC440
  blas_ಸ್ಮರಣೆ_free(ಬಫರ್);
#ಆದರೆ_ಅಂತ್ಯ

  *Info = info;

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;
}
