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

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ SSCAL
  FLOAT *alpha=ALPHA;
#ਨਹੀਂਤਾਂ
  FLOAT alpha[2] = {ALPHA[0], ZERO};
#ਜੇ_ਅੰਤ

#ਨਹੀਂਤਾਂ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ SSCAL
ਖਾਲੀ CNAME(blasint n, ਖਾਲੀ *VALPHA, ਖਾਲੀ *vx, blasint incx){

  FLOAT *x = (FLOAT*) vx;
  FLOAT *alpha=(FLOAT*)VALPHA;
#ਨਹੀਂਤਾਂ
ਖਾਲੀ CNAME(blasint n, FLOAT alpha_r, ਖਾਲੀ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;

  FLOAT alpha[2] = {alpha_r, ZERO};
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ ਬਹੁਲਕ;
  ਪੂਰਨ nthreads;
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
  PRINT_DEBUG_NAME;
#ਨਹੀਂਤਾਂ
  PRINT_DEBUG_CNAME;
#ਜੇ_ਅੰਤ

  ਜੇ (incx <= 0 || n <= 0) ਵਾਪਸ;

  ਜੇ ((alpha[0] == ONE) && (alpha[1] == ZERO)) ਵਾਪਸ;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਜੇ ( n <= 1048576 )
	nthreads = 1;
  ਨਹੀਂਤਾਂ
	nthreads = ਅੰਕ_cpu_avail(1);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  SCAL_K(n, 0, 0, alpha[0], alpha[1], x, incx, ਨਿਰਾ, 0, ਨਿਰਾ, 1);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
    ਬਹੁਲਕ  =  BLAS_DOUBLE | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
    ਬਹੁਲਕ  =  BLAS_SINGLE | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ

    blas_level1_thread(ਬਹੁਲਕ, n, 0, 0,  alpha, x, incx, ਨਿਰਾ, 0, ਨਿਰਾ, 1, (ਪੂਰਨ (*)(ਖਾਲੀ))SCAL_K, nthreads);

  }
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(4, n, n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;

}
