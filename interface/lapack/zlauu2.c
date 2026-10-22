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
#تعریف_کن ERROR_NAME "QLAUU2"
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZLAUU2"
#وإلا
#تعریف_کن ERROR_NAME "CLAUU2"
#پایان_اگر

ساكن blasint (*lauu2[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {

#اگر_تعریف_شده XDOUBLE
  xlauu2_U, xlauu2_L,
#وإلاإذا defined(DOUBLE)
  zlauu2_U, zlauu2_L,
#وإلا
  clauu2_U, clauu2_L,
#پایان_اگر
  };

OPENBLAS_EXPORT
فراغ NAME(حرف *UPLO, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t آرگومانها;

  blasint uplo_arg = *UPLO;
  blasint uplo;
  blasint info;
  FLOAT *میانگیر;
#اگر_تعریف_شده PPC440
  خارجي
#پایان_اگر
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  آرگومانها.n    = *N;
  آرگومانها.a    = (فراغ *)a;
  آرگومانها.lda  = *ldA;

  TOUPPER(uplo_arg);

  uplo = -1;
  إذا (uplo_arg == 'U') uplo = 0;
  إذا (uplo_arg == 'L') uplo = 1;

  info  = 0;
  إذا (آرگومانها.lda < MAX(1,آرگومانها.n)) info = 4;
  إذا (آرگومانها.n   < 0)             info = 2;
  إذا (uplo     < 0)             info = 1;
  إذا (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME) - 1);
    *Info = - info;
    أرجع;
  }

  *Info = 0;

  إذا (آرگومانها.n <= 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

#اگر_تعریف_نشده PPC440
  میانگیر = (FLOAT *)blas_حافظه_alloc(1);

  sa = (FLOAT *)((BLASLONG)میانگیر + GEMM_جابه_جایی_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * اندازه + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_جابه_جایی_B);
#پایان_اگر

  info = (lauu2[uplo])(&آرگومانها, عدم, عدم, sa, sb, 0);

  *Info = info;

#اگر_تعریف_نشده PPC440
  blas_حافظه_free(میانگیر);
#پایان_اگر

  FUNCTION_PROFILE_انتها(1, .5 * آرگومانها.n * آرگومانها.n,
		       2. * آرگومانها.n * (1./3. + آرگومانها.n * ( 1./2. + آرگومانها.n * 1./6.))
		       +  6. * 1./6. * آرگومانها.n * (آرگومانها.n * آرگومانها.n - 1));

  IDEBUG_انتها;

  أرجع;
}
