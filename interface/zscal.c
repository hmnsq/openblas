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

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;

#ifndef SSCAL
  FLOAT *alpha=ALPHA;
#അല്ലെങ്കിൽ
  FLOAT alpha[2] = {ALPHA[0], ZERO};
#endif

#അല്ലെങ്കിൽ

#ifndef SSCAL
ശൂന്യം CNAME(blasint n, ശൂന്യം *VALPHA, ശൂന്യം *vx, blasint incx){

  FLOAT *x = (FLOAT*) vx;
  FLOAT *alpha=(FLOAT*)VALPHA;
#അല്ലെങ്കിൽ
ശൂന്യം CNAME(blasint n, FLOAT alpha_r, ശൂന്യം *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;

  FLOAT alpha[2] = {alpha_r, ZERO};
#endif
#endif

#ifdef SMP
  പൂർണ്ണസംഖ്യ mode;
  പൂർണ്ണസംഖ്യ nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#അല്ലെങ്കിൽ
  PRINT_DEBUG_CNAME;
#endif

  എങ്കിൽ (incx <= 0 || n <= 0) തിരികെ;

  എങ്കിൽ ((alpha[0] == ONE) && (alpha[1] == ZERO)) തിരികെ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifdef SMP
  എങ്കിൽ ( n <= 1048576 )
	nthreads = 1;
  അല്ലെങ്കിൽ
	nthreads = അക്കം_cpu_avail(1);

  എങ്കിൽ (nthreads == 1) {
#endif

  SCAL_K(n, 0, 0, alpha[0], alpha[1], x, incx, ശൂന്യ, 0, ശൂന്യ, 1);

#ifdef SMP
  } അല്ലെങ്കിൽ {
#ifdef DOUBLE
    mode  =  BLAS_DOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
    mode  =  BLAS_SINGLE | BLAS_COMPLEX;
#endif

    blas_level1_thread(mode, n, 0, 0,  alpha, x, incx, ശൂന്യ, 0, ശൂന്യ, 1, (പൂർണ്ണസംഖ്യ (*)(ശൂന്യം))SCAL_K, nthreads);

  }
#endif

  FUNCTION_PROFILE_END(4, n, n);

  IDEBUG_END;

  തിരികെ;

}
