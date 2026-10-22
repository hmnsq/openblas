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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#إذا defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#وإلا
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#وصف_رد SMP
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#وإلا

فراغ CNAME(blasint n, فراغ *vx, blasint incx, فراغ *vy, blasint incy){
FLOAT *x = (FLOAT*)vx;
FLOAT *y = (FLOAT*)vy;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
  صحيح گهڻائي;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  صحيح nthreads;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جيڪڏهن_پورو

  إذا (n <= 0) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#جيڪڏهن_وصف_ٿيل SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  إذا (incx == 0 || incy == 0 || n < 1048576 * GEMM_MULTITHREAD_حد / حجم(FLOAT))
	  nthreads = 1;
  وإلا
	  nthreads = انگ_cpu_avail(1);

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

  مٽا_سٽا_K(n, 0, 0, ZERO, ZERO, x, incx, y, incy, عدم, 0);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

#جيڪڏهن_وصف_ٿيل XDOUBLE
    گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
    گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
    گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو

    blas_level1_thread(گهڻائي, n, 0, 0, dummyalpha,
		       x, incx, y, incy, عدم, 0, (فراغ *)مٽا_سٽا_K, nthreads);

  }
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(2, 2 * n, 0);

  IDEBUG_پڇاڙي;

  أرجع;

}
