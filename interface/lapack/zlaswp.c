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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

സ്ഥിതം പൂർണ്ണസംഖ്യ (*laswp[])(BLASLONG, BLASLONG, BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, blasint *, BLASLONG) = {
#ifdef XDOUBLE
  xlaswp_plus, xlaswp_minus,
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  zlaswp_plus, zlaswp_minus,
#അല്ലെങ്കിൽ
  claswp_plus, claswp_minus,
#endif
};

OPENBLAS_EXPORT
ശൂന്യം NAME(blasint *N, FLOAT *a, blasint *LDA, blasint *K1, blasint *K2, blasint *ipiv, blasint *INCX){

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint k1   = *K1;
  blasint k2   = *K2;
  blasint incx = *INCX;
  പൂർണ്ണസംഖ്യ flag;

#ifdef SMP
  പൂർണ്ണസംഖ്യ mode;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  പൂർണ്ണസംഖ്യ nthreads;
#endif

  PRINT_DEBUG_NAME;

  എങ്കിൽ (incx == 0 || n <= 0) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  flag = (incx < 0);

#ifdef SMP
  എങ്കിൽ (n < 64) {
    nthreads = 1;
  } അല്ലെങ്കിൽ {
    nthreads = അക്കം_cpu_avail(2);
  }

  എങ്കിൽ (nthreads == 1) {
#endif

  (laswp[flag])(n, k1, k2, ZERO, ZERO, a, lda, ശൂന്യ, 0, ipiv, incx);

#ifdef SMP
  } അല്ലെങ്കിൽ {

#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif

  blas_level1_thread(mode, n, k1, k2, dummyalpha, a, lda, ശൂന്യ, 0, ipiv, incx, (പൂർണ്ണസംഖ്യ(*)(ശൂന്യം))laswp[flag], nthreads);
  }
#endif

  FUNCTION_PROFILE_END(COMPSIZE, n * (k2 - k1), 0);

  IDEBUG_END;

  തിരികെ;
}
