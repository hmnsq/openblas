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

#ضمن <stdio.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

#إذا defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#وإلا
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#لغو_تعریف SMP
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#وإلا

فراغ CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#پایان_اگر

#اگر_تعریف_شده SMP
  صحيح نما, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#پایان_اگر

#اگر_تعریف_نشده CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#پایان_اگر

  إذا (n <= 0) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0) x -= (n - 1) * incx;
  إذا (incy < 0) y -= (n - 1) * incy;

#اگر_تعریف_شده SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  إذا (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_آستانه / حجم(FLOAT))
    nthreads = 1;
  وإلا
    nthreads = رقم_cpu_avail(1);

  إذا (nthreads == 1) {
#پایان_اگر

    جابجا_K(n, 0, 0, ZERO, x, incx, y, incy, عدم, 0);

#اگر_تعریف_شده SMP
  } وإلا {

#اگر_تعریف_شده XDOUBLE
    نما  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
    نما  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
    نما  =  BLAS_SINGLE  | BLAS_REAL;
#پایان_اگر

    blas_level1_thread(نما, n, 0, 0, dummyalpha,
		       x, incx, y, incy, عدم, 0, (فراغ *)جابجا_K, nthreads);
  }

#پایان_اگر

  FUNCTION_PROFILE_انتها(1, 2 * n, 0);

  IDEBUG_انتها;

  أرجع;

}
