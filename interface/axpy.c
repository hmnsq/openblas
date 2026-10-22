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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای
#إذا  defined(Z13)
#تعریف_کړه MULTI_THREAD_اصغري  200000
#وإلا
#تعریف_کړه MULTI_THREAD_اصغري  10000
#د_که_پای

OPENBLAS_EXPORT
#که_نه_تعریف_شوی CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT alpha = *ALPHA;

#وإلا

فراغ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#د_که_پای

#که_تعریف_شوی SMP
  صحيح ډېرپېښ, nthreads;
#د_که_پای

#که_نه_تعریف_شوی CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#د_که_پای

  إذا (n <= 0) أرجع;

  إذا (alpha == ZERO) أرجع;

  إذا (incx == 0 && incy == 0) {
    *y += n * alpha *(*x);
    أرجع;
  }
    
  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  إذا (incx < 0) x -= (n - 1) * incx;
  إذا (incy < 0) y -= (n - 1) * incy;

#که_تعریف_شوی SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  إذا (incx == 0 || incy == 0 || n <= MULTI_THREAD_اصغري)
	  nthreads = 1;
  وإلا
	  nthreads = عدد_cpu_avail(1);

  إذا (nthreads == 1) {
#د_که_پای

  AXPYU_K(n, 0, 0, alpha, x, incx, y, incy, عدم, 0);

#که_تعریف_شوی SMP
  } وإلا {

#که_تعریف_شوی XDOUBLE
    ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
    ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
    ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای

    blas_level1_thread(ډېرپېښ, n, 0, 0, &alpha,
		       x, incx, y, incy, عدم, 0,  (صحيح (*)(فراغ))AXPYU_K, nthreads);

  }
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(1, 2 * n, 2 * n);

  IDEBUG_پای_ټکی;

  أرجع;

}
