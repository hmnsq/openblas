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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif
#જો  defined(Z13)
#define MULTI_THREAD_MINIMAL  200000
#નહીંતર
#define MULTI_THREAD_MINIMAL  10000
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#નહીંતર
#ifdef COMPLEX
શૂન્ય CNAME(blasint n, શૂન્ય *VALPHA, શૂન્ય *vx, blasint incx, શૂન્ય *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#નહીંતર
શૂન્ય CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#endif
#endif

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#ifdef SMP
  પૂર્ણાંક mode, nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#નહીંતર
  PRINT_DEBUG_CNAME;
#endif

  જો (n <= 0) પાછા;

  જો ((alpha_r == ZERO) && (alpha_i == ZERO)) પાછા;

  જો (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  પાછા;
  }
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (n - 1) * incx * 2;
  જો (incy < 0) y -= (n - 1) * incy * 2;

#ifdef SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  જો (incx == 0 || incy == 0 || n <= MULTI_THREAD_MINIMAL)
	  nthreads = 1;
  નહીંતર
	  nthreads = અંક_cpu_avail(1);

  જો (nthreads == 1) {
#endif

#ifndef CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, ખાલી, 0);
#નહીંતર
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, ખાલી, 0);
#endif

#ifdef SMP
  } નહીંતર {

#ifdef XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
    mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif

    blas_level1_thread(mode, n, 0, 0, ALPHA, x, incx, y, incy, ખાલી, 0,
#ifndef CONJ
                       (પૂર્ણાંક (*)(શૂન્ય))AXPYU_K,
#નહીંતર
                       (પૂર્ણાંક (*)(શૂન્ય))AXPYC_K,
#endif
		       nthreads);
  }
#endif

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  IDEBUG_END;

  પાછા;

}
