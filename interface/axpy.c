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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर
#अगर  defined(Z13)
#घोषणा MULTI_THREAD_MINIMAL  200000
#अन्यथा
#घोषणा MULTI_THREAD_MINIMAL  10000
#पूर्ण_अगर

OPENBLAS_EXPORT
#अगर_अघोषित CBLAS

व्योम NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT alpha = *ALPHA;

#अन्यथा

व्योम CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#पूर्ण_अगर

#अगर_घोषित SMP
  पूर्णांक mode, nthreads;
#पूर्ण_अगर

#अगर_अघोषित CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#पूर्ण_अगर

  अगर (n <= 0) वापस;

  अगर (alpha == ZERO) वापस;

  अगर (incx == 0 && incy == 0) {
    *y += n * alpha *(*x);
    वापस;
  }
    
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (n - 1) * incx;
  अगर (incy < 0) y -= (n - 1) * incy;

#अगर_घोषित SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  अगर (incx == 0 || incy == 0 || n <= MULTI_THREAD_MINIMAL)
	  nthreads = 1;
  अन्यथा
	  nthreads = अंक_cpu_avail(1);

  अगर (nthreads == 1) {
#पूर्ण_अगर

  AXPYU_K(n, 0, 0, alpha, x, incx, y, incy, शून्य, 0);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_घोषित XDOUBLE
    mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
    mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
    mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर

    blas_level1_thread(mode, n, 0, 0, &alpha,
		       x, incx, y, incy, शून्य, 0,  (पूर्णांक (*)(व्योम))AXPYU_K, nthreads);

  }
#पूर्ण_अगर

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);

  IDEBUG_END;

  वापस;

}
