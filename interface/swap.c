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
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ಯದಿ defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#ಅನ್ಯಥಾ
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#undef SMP
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#endif

#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ mode, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ಅನ್ಯಥಾ
  PRINT_DEBUG_CNAME;
#endif

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (n - 1) * incx;
  ಯದಿ (incy < 0) y -= (n - 1) * incy;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  ಯದಿ (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_THRESHOLD / ಪರಿಮಾಣ(FLOAT))
    nthreads = 1;
  ಅನ್ಯಥಾ
    nthreads = ಅಂಕೆ_cpu_avail(1);

  ಯದಿ (nthreads == 1) {
#endif

    SWAP_K(n, 0, 0, ZERO, x, incx, y, incy, ರಿಕ್ತ, 0);

#ifdef SMP
  } ಅನ್ಯಥಾ {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
    mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif

    blas_level1_thread(mode, n, 0, 0, dummyalpha,
		       x, incx, y, incy, ರಿಕ್ತ, 0, (ಶೂನ್ಯ *)SWAP_K, nthreads);
  }

#endif

  FUNCTION_PROFILE_END(1, 2 * n, 0);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;

}
