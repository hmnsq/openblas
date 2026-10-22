/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project.                              */
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

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#अन्यथा

शून्य CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#यदि_अन्तः

#यदि_परिभाषितम् SMP
  पूर्णाङ्क बहुलकम्, nthreads;
#यदि_अन्तः

#यदि_अपरिभाषितम् CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#यदि_अन्तः

  यदि (incx <= 0 || n <= 0) प्रत्यागम;

#यदि_परिभाषितम् BGEMM
  भिन्न alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#अन्यथा
  FLOAT alpha_float = alpha;
#यदि_अन्तः

  यदि (alpha_float == ONE) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();


#यदि_परिभाषितम् SMP
  यदि (n <= 1048576 )
	nthreads = 1;
  अन्यथा
	nthreads = अङ्कः_cpu_avail(1);


  यदि (nthreads == 1) {
#यदि_अन्तः

  SCAL_K(n, 0, 0, alpha, x, incx, रिक्त, 0, रिक्त, 1);

#यदि_परिभाषितम् SMP
  } अन्यथा {

#यदि_परिभाषितम् DOUBLE
    बहुलकम्  =  BLAS_DOUBLE | BLAS_REAL;
#अन्यथा
    बहुलकम्  =  BLAS_SINGLE | BLAS_REAL;
#यदि_अन्तः

    blas_level1_thread(बहुलकम्, n, 0, 0,
#यदि_अपरिभाषितम् CBLAS
		       ALPHA,
#अन्यथा
		       &alpha,
#यदि_अन्तः
		       x, incx, रिक्त, 0, रिक्त, 1,  (पूर्णाङ्क (*)(शून्य))SCAL_K, nthreads);

  }
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(1, n, n);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
