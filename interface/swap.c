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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#अन्यथा
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#अपरिभाषय SMP
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#अन्यथा

शून्य CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#यदि_अन्तः

#यदि_परिभाषितम् SMP
  पूर्णाङ्क बहुलकम्, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#यदि_अन्तः

#यदि_अपरिभाषितम् CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#यदि_अन्तः

  यदि (n <= 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (n - 1) * incx;
  यदि (incy < 0) y -= (n - 1) * incy;

#यदि_परिभाषितम् SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  यदि (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_देहली / परिमाण(FLOAT))
    nthreads = 1;
  अन्यथा
    nthreads = अङ्कः_cpu_avail(1);

  यदि (nthreads == 1) {
#यदि_अन्तः

    SWAP_K(n, 0, 0, ZERO, x, incx, y, incy, रिक्त, 0);

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_परिभाषितम् XDOUBLE
    बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
    बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
    बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः

    blas_level1_thread(बहुलकम्, n, 0, 0, dummyalpha,
		       x, incx, y, incy, रिक्त, 0, (शून्य *)SWAP_K, nthreads);
  }

#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(1, 2 * n, 0);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
