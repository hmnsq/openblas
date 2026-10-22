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
#إذا  defined(Z13)
#وصف_ڪريو MULTI_THREAD_گهٽتر  200000
#وإلا
#وصف_ڪريو MULTI_THREAD_گهٽتر  10000
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#وإلا
#جيڪڏهن_وصف_ٿيل COMPLEX
فراغ CNAME(blasint n, فراغ *VALPHA, فراغ *vx, blasint incx, فراغ *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#وإلا
فراغ CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#جيڪڏهن_وصف_ٿيل SMP
  صحيح گهڻائي, nthreads;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جيڪڏهن_پورو

  إذا (n <= 0) أرجع;

  إذا ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  إذا (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  أرجع;
  }
  
  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#جيڪڏهن_وصف_ٿيل SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  إذا (incx == 0 || incy == 0 || n <= MULTI_THREAD_گهٽتر)
	  nthreads = 1;
  وإلا
	  nthreads = انگ_cpu_avail(1);

  إذا (nthreads == 1) {
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, عدم, 0);
#وإلا
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, عدم, 0);
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

#جيڪڏهن_وصف_ٿيل XDOUBLE
    گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
    گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
    گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو

    blas_level1_thread(گهڻائي, n, 0, 0, ALPHA, x, incx, y, incy, عدم, 0,
#جيڪڏهن_وصف_نه_ٿيل CONJ
                       (صحيح (*)(فراغ))AXPYU_K,
#وإلا
                       (صحيح (*)(فراغ))AXPYC_K,
#جيڪڏهن_پورو
		       nthreads);
  }
#جيڪڏهن_پورو

  FUNCTION_PROFILE_پڇاڙي(4, 2 * n, 2 * n);

  IDEBUG_پڇاڙي;

  أرجع;

}
