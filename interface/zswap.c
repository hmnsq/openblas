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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#nếu defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#không_thì
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#undef SMP
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#không_thì

void CNAME(blasint n, void *vx, blasint incx, void *vy, blasint incy){
FLOAT *x = (FLOAT*)vx;
FLOAT *y = (FLOAT*)vy;
#endif

#ifdef SMP
  số_nguyên mode;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  số_nguyên nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#không_thì
  PRINT_DEBUG_CNAME;
#endif

  nếu (n <= 0) trả_về;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  nếu (incx < 0) x -= (n - 1) * incx * 2;
  nếu (incy < 0) y -= (n - 1) * incy * 2;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  nếu (incx == 0 || incy == 0 || n < 1048576 * GEMM_MULTITHREAD_THRESHOLD / sizeof(FLOAT))
	  nthreads = 1;
  không_thì
	  nthreads = num_cpu_avail(1);

  nếu (nthreads == 1) {
#endif

  SWAP_K(n, 0, 0, ZERO, ZERO, x, incx, y, incy, NULL, 0);

#ifdef SMP
  } không_thì {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#không_thì
    mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif

    blas_level1_thread(mode, n, 0, 0, dummyalpha,
		       x, incx, y, incy, NULL, 0, (void *)SWAP_K, nthreads);

  }
#endif

  FUNCTION_PROFILE_END(2, 2 * n, 0);

  IDEBUG_END;

  trả_về;

}
