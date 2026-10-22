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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#सचे defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#अञ्ञथा
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#undef SMP
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#अञ्ञथा

शून्य CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#endif

#ifdef SMP
  पूर्णाङ्क mode, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#अञ्ञथा
  PRINT_DEBUG_CNAME;
#endif

  सचे (n <= 0) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0) x -= (n - 1) * incx;
  सचे (incy < 0) y -= (n - 1) * incy;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  सचे (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_THRESHOLD / परिमाण(FLOAT))
    nthreads = 1;
  अञ्ञथा
    nthreads = अङ्कः_cpu_avail(1);

  सचे (nthreads == 1) {
#endif

    SWAP_K(n, 0, 0, ZERO, x, incx, y, incy, रिक्त, 0);

#ifdef SMP
  } अञ्ञथा {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अञ्ञथा
    mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif

    blas_level1_thread(mode, n, 0, 0, dummyalpha,
		       x, incx, y, incy, रिक्त, 0, (शून्य *)SWAP_K, nthreads);
  }

#endif

  FUNCTION_PROFILE_END(1, 2 * n, 0);

  IDEBUG_END;

  पटिनिवत्त;

}
