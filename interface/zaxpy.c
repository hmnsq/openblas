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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت
#إذا  defined(Z13)
#تعریف_کرو MULTI_THREAD_اقل  200000
#وإلا
#تعریف_کرو MULTI_THREAD_اقل  10000
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#وإلا
#جے_تعریف_شدہ COMPLEX
فراغ CNAME(blasint n, فراغ *VALPHA, فراغ *vx, blasint incx, فراغ *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#وإلا
فراغ CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#جے_انت
#جے_انت

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#جے_تعریف_شدہ SMP
  صحيح منوال, nthreads;
#جے_انت

#جے_بنا_تعریف CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جے_انت

  إذا (n <= 0) أرجع;

  إذا ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  إذا (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  أرجع;
  }
  
  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  إذا (incx < 0) x -= (n - 1) * incx * 2;
  إذا (incy < 0) y -= (n - 1) * incy * 2;

#جے_تعریف_شدہ SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  إذا (incx == 0 || incy == 0 || n <= MULTI_THREAD_اقل)
	  nthreads = 1;
  وإلا
	  nthreads = عدد_cpu_avail(1);

  إذا (nthreads == 1) {
#جے_انت

#جے_بنا_تعریف CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, عدم, 0);
#وإلا
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, عدم, 0);
#جے_انت

#جے_تعریف_شدہ SMP
  } وإلا {

#جے_تعریف_شدہ XDOUBLE
    منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
    منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
    منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت

    blas_level1_thread(منوال, n, 0, 0, ALPHA, x, incx, y, incy, عدم, 0,
#جے_بنا_تعریف CONJ
                       (صحيح (*)(فراغ))AXPYU_K,
#وإلا
                       (صحيح (*)(فراغ))AXPYC_K,
#جے_انت
		       nthreads);
  }
#جے_انت

  FUNCTION_PROFILE_سرا(4, 2 * n, 2 * n);

  IDEBUG_سرا;

  أرجع;

}
