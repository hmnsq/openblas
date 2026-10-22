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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः
#यदि  defined(Z13)
#परिभाषय MULTI_THREAD_अल्पिष्ठम्  200000
#अन्यथा
#परिभाषय MULTI_THREAD_अल्पिष्ठम्  10000
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#अन्यथा
#यदि_परिभाषितम् COMPLEX
शून्य CNAME(blasint n, शून्य *VALPHA, शून्य *vx, blasint incx, शून्य *vy, blasint incy){
FLOAT *ALPHA = (FLOAT*) VALPHA;
FLOAT *x = (FLOAT*) vx;
FLOAT *y = (FLOAT*) vy;
#अन्यथा
शून्य CNAME(blasint n, FLOAT *ALPHA, FLOAT *x, blasint incx, FLOAT *y, blasint incy){
#यदि_अन्तः
#यदि_अन्तः

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);

#यदि_परिभाषितम् SMP
  पूर्णाङ्क बहुलकम्, nthreads;
#यदि_अन्तः

#यदि_अपरिभाषितम् CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#यदि_अन्तः

  यदि (n <= 0) प्रत्यागम;

  यदि ((alpha_r == ZERO) && (alpha_i == ZERO)) प्रत्यागम;

  यदि (incx == 0 && incy == 0) {
  *y += n * (alpha_r * (*x) - alpha_i* (*(x+1)) );
  *(y+1) += n * (alpha_i * (*x) + alpha_r * (*(x +1)) );
  प्रत्यागम;
  }
  
  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (n - 1) * incx * 2;
  यदि (incy < 0) y -= (n - 1) * incy * 2;

#यदि_परिभाषितम् SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  //
  //Temporarily work-around the low performance issue with small input size &
  //multithreads.
  यदि (incx == 0 || incy == 0 || n <= MULTI_THREAD_अल्पिष्ठम्)
	  nthreads = 1;
  अन्यथा
	  nthreads = अङ्कः_cpu_avail(1);

  यदि (nthreads == 1) {
#यदि_अन्तः

#यदि_अपरिभाषितम् CONJ
    AXPYU_K (n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, रिक्त, 0);
#अन्यथा
    AXPYC_K(n, 0, 0, alpha_r, alpha_i, x, incx, y, incy, रिक्त, 0);
#यदि_अन्तः

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_परिभाषितम् XDOUBLE
    बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
    बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
    बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः

    blas_level1_thread(बहुलकम्, n, 0, 0, ALPHA, x, incx, y, incy, रिक्त, 0,
#यदि_अपरिभाषितम् CONJ
                       (पूर्णाङ्क (*)(शून्य))AXPYU_K,
#अन्यथा
                       (पूर्णाङ्क (*)(शून्य))AXPYC_K,
#यदि_अन्तः
		       nthreads);
  }
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(4, 2 * n, 2 * n);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
