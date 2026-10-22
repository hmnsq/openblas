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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#अन्यथा
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#अघोषित SMP
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#अन्यथा

व्योम CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#पूर्ण_अगर

#अगर_घोषित SMP
  पूर्णांक mode, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#पूर्ण_अगर

#अगर_अघोषित CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#पूर्ण_अगर

  अगर (n <= 0) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (n - 1) * incx;
  अगर (incy < 0) y -= (n - 1) * incy;

#अगर_घोषित SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  अगर (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_THRESHOLD / माप(FLOAT))
    nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(1);

  अगर (nthreads == 1) {
#पूर्ण_अगर

    SWAP_K(n, 0, 0, ZERO, x, incx, y, incy, शून्य, 0);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_घोषित XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
    mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर

    blas_level1_thread(mode, n, 0, 0, dummyalpha,
		       x, incx, y, incy, शून्य, 0, (व्योम *)SWAP_K, nthreads);
  }

#पूर्ण_अगर

  FUNCTION_PROFILE_END(1, 2 * n, 0);

  IDEBUG_END;

  वापस;

}
