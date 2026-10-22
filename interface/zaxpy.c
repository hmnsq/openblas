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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif
#যদি  defined(Z13)
#define MULTI_THREAD_MINIMAL  200000
#অন্যথা
#define MULTI_THREAD_MINIMAL  10000
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#অন্যথা
#ifdef COMPLEX
শূন্য CNAME(blasint n, শূন্য *VALPHA, শূন্য *vx, blasint incx, শূন্য *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#অন্যথা
শূন্য CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#endif
#endif

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#ifdef SMP
  পূর্ণাঙ্ক mode, nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#অন্যথা
  PRINT_DEBUG_CNAME;
#endif

  যদি (n <= 0) ঘূৰাই;

  যদি ((alpha_r == ZERO) && (alpha_i == ZERO)) ঘূৰাই;

  যদি (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  ঘূৰাই;
  }
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (n - 1) * incx * 2;
  যদি (incy < 0) y -= (n - 1) * incy * 2;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  যদি (incx == 0 || incy == 0 || n <= MULTI_THREAD_MINIMAL)
	  nthreads = 1;
  অন্যথা
	  nthreads = অঙ্ক_cpu_avail(1);

  যদি (nthreads == 1) {
#endif

#ifndef CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, রিক্ত, 0);
#অন্যথা
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, রিক্ত, 0);
#endif

#ifdef SMP
  } অন্যথা {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
    mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif

    blas_level1_thread(mode, n, 0, 0, ALPHA, x, incx, y, incy, রিক্ত, 0,
#ifndef CONJ
                       (পূর্ণাঙ্ক (*)(শূন্য))AXPYU_K,
#অন্যথা
                       (পূর্ণাঙ্ক (*)(শূন্য))AXPYC_K,
#endif
		       nthreads);
  }
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  ঘূৰাই;

}
