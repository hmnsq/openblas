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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XTRTI2"
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZTRTI2"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CTRTI2"
#ਜੇ_ਅੰਤ

ਸਥਿਰ blasint (*trti2[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  xtrti2_UU, xtrti2_UN, xtrti2_LU, xtrti2_LN,
#ਜਾਂਜੇ defined(DOUBLE)
  ztrti2_UU, ztrti2_UN, ztrti2_LU, ztrti2_LN,
#ਨਹੀਂਤਾਂ
  ctrti2_UU, ctrti2_UN, ctrti2_LU, ctrti2_LN,
#ਜੇ_ਅੰਤ
  };

OPENBLAS_EXPORT
ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *DIAG, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t ਦਲੀਲਾਂ;

  blasint uplo_arg = *UPLO;
  blasint diag_arg = *DIAG;
  blasint uplo, diag;
  blasint info;
  FLOAT *ਬਫਰ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ PPC440
  ਬਾਹਰੀ
#ਜੇ_ਅੰਤ
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  ਦਲੀਲਾਂ.n    = *N;
  ਦਲੀਲਾਂ.a    = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.lda  = *ldA;

  TOUPPER(uplo_arg);
  TOUPPER(diag_arg);

  uplo = -1;
  ਜੇ (uplo_arg == 'U') uplo = 0;
  ਜੇ (uplo_arg == 'L') uplo = 1;
  diag = -1;
  ਜੇ (diag_arg == 'U') diag = 0;
  ਜੇ (diag_arg == 'N') diag = 1;

  info  = 0;
  ਜੇ (ਦਲੀਲਾਂ.lda  < MAX(1,ਦਲੀਲਾਂ.n)) info = 5;
  ਜੇ (ਦਲੀਲਾਂ.n    < 0)             info = 3;
  ਜੇ (diag < 0)                  info = 2;
  ਜੇ (uplo < 0)                  info = 1;
  ਜੇ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME) - 1);
    *Info = - info;
    ਵਾਪਸ;
  }

  *Info = 0;

  ਜੇ (ਦਲੀਲਾਂ.n <= 0) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ PPC440
  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(1);

  sa = (FLOAT *)((BLASLONG)ਬਫਰ + GEMM_ਵਿੱਥ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);
#ਜੇ_ਅੰਤ

  info = (trti2[(uplo << 1) | diag])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

  *Info = info;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ PPC440
  blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(1, .5 * ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.n,
		          2. * ਦਲੀਲਾਂ.n * (1./3. + ਦਲੀਲਾਂ.n * ( 1./2. + ਦਲੀਲਾਂ.n * 1./6.))
		       +  6. * ਦਲੀਲਾਂ.n * (1./3. + ਦਲੀਲਾਂ.n * (-1./2. + ਦਲੀਲਾਂ.n * 1./6.)));

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
