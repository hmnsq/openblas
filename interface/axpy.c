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
#ਜੇ  defined(Z13)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MULTI_THREAD_ਨਿਊਨਤਮ  200000
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MULTI_THREAD_ਨਿਊਨਤਮ  10000
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT alpha = *ALPHA;

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ ਬਹੁਲਕ, nthreads;
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
  PRINT_DEBUG_NAME;
#ਨਹੀਂਤਾਂ
  PRINT_DEBUG_CNAME;
#ਜੇ_ਅੰਤ

  ਜੇ (n <= 0) ਵਾਪਸ;

  ਜੇ (alpha == ZERO) ਵਾਪਸ;

  ਜੇ (incx == 0 && incy == 0) {
    *y += n * alpha *(*x);
    ਵਾਪਸ;
  }
    
  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਜੇ (incx < 0) x -= (n - 1) * incx;
  ਜੇ (incy < 0) y -= (n - 1) * incy;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  ਜੇ (incx == 0 || incy == 0 || n <= MULTI_THREAD_ਨਿਊਨਤਮ)
	  nthreads = 1;
  ਨਹੀਂਤਾਂ
	  nthreads = ਅੰਕ_cpu_avail(1);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ

  AXPYU_K(n, 0, 0, alpha, x, incx, y, incy, ਨਿਰਾ, 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
    ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
    ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
    ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ

    blas_level1_thread(ਬਹੁਲਕ, n, 0, 0, &alpha,
		       x, incx, y, incy, ਨਿਰਾ, 0,  (ਪੂਰਨ (*)(ਖਾਲੀ))AXPYU_K, nthreads);

  }
#ਜੇ_ਅੰਤ

  FUNCTION_PROFILE_ਸਿਰਾ(1, 2 * n, 2 * n);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;

}
