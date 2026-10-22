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

#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم
#اگر  defined(Z13)
#تعریف_کرو MULTI_THREAD_اقل  200000
#ورنہ
#تعریف_کرو MULTI_THREAD_اقل  10000
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#ورنہ
#اگر_تعریف_شدہ COMPLEX
خالی CNAME(blasint n, خالی *VALPHA, خالی *vx, blasint incx, خالی *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#ورنہ
خالی CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#اگر_ختم
#اگر_ختم

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#اگر_تعریف_شدہ SMP
  صحیح منوال, nthreads;
#اگر_ختم

#اگر_غیر_تعریف_شدہ CBLAS
  PRINT_DEBUG_NAME;
#ورنہ
  PRINT_DEBUG_CNAME;
#اگر_ختم

  اگر (n <= 0) واپس;

  اگر ((alpha_r == ZERO) && (alpha_i == ZERO)) واپس;

  اگر (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  واپس;
  }
  
  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0) x -= (n - 1) * incx * 2;
  اگر (incy < 0) y -= (n - 1) * incy * 2;

#اگر_تعریف_شدہ SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  اگر (incx == 0 || incy == 0 || n <= MULTI_THREAD_اقل)
	  nthreads = 1;
  ورنہ
	  nthreads = عدد_cpu_avail(1);

  اگر (nthreads == 1) {
#اگر_ختم

#اگر_غیر_تعریف_شدہ CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, تہی, 0);
#ورنہ
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, تہی, 0);
#اگر_ختم

#اگر_تعریف_شدہ SMP
  } ورنہ {

#اگر_تعریف_شدہ XDOUBLE
    منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
    منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
    منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم

    blas_level1_thread(منوال, n, 0, 0, ALPHA, x, incx, y, incy, تہی, 0,
#اگر_غیر_تعریف_شدہ CONJ
                       (صحیح (*)(خالی))AXPYU_K,
#ورنہ
                       (صحیح (*)(خالی))AXPYC_K,
#اگر_ختم
		       nthreads);
  }
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, 2 * n, 2 * n);

  IDEBUG_انتہا;

  واپس;

}
