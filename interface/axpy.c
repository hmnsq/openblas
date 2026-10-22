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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif
#ଯଦି  defined(Z13)
#define MULTI_THREAD_MINIMAL  200000
#ଅନ୍ୟଥା
#define MULTI_THREAD_MINIMAL  10000
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

ଶୂନ୍ୟ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT alpha = *ALPHA;

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#endif

#ifdef SMP
  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode, nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ଅନ୍ୟଥା
  PRINT_DEBUG_CNAME;
#endif

  ଯଦି (n <= 0) ଫେରନ୍ତୁ;

  ଯଦି (alpha == ZERO) ଫେରନ୍ତୁ;

  ଯଦି (incx == 0 && incy == 0) {
    *y += n * alpha *(*x);
    ଫେରନ୍ତୁ;
  }
    
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ଯଦି (incx < 0) x -= (n - 1) * incx;
  ଯଦି (incy < 0) y -= (n - 1) * incy;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  ଯଦି (incx == 0 || incy == 0 || n <= MULTI_THREAD_MINIMAL)
	  nthreads = 1;
  ଅନ୍ୟଥା
	  nthreads = ଅଙ୍କ_cpu_avail(1);

  ଯଦି (nthreads == 1) {
#endif

  AXPYU_K(n, 0, 0, alpha, x, incx, y, incy, ଖାଲି, 0);

#ifdef SMP
  } ଅନ୍ୟଥା {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
    mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif

    blas_level1_thread(mode, n, 0, 0, &alpha,
		       x, incx, y, incy, ଖାଲି, 0,  (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(ଶୂନ୍ୟ))AXPYU_K, nthreads);

  }
#endif

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);

  IDEBUG_END;

  ଫେରନ୍ତୁ;

}
