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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#अन्यथा

व्योम CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#पूर्ण_अगर

#अगर_घोषित SMP
  पूर्णांक mode, nthreads;
#पूर्ण_अगर

#अगर_अघोषित CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#पूर्ण_अगर

  अगर (incx <= 0 || n <= 0) वापस;

#अगर_घोषित BGEMM
  भग्न alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#अन्यथा
  FLOAT alpha_float = alpha;
#पूर्ण_अगर

  अगर (alpha_float == ONE) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();


#अगर_घोषित SMP
  अगर (n <= 1048576 )
	nthreads = 1;
  अन्यथा
	nthreads = अंक_cpu_avail(1);


  अगर (nthreads == 1) {
#पूर्ण_अगर

  SCAL_K(n, 0, 0, alpha, x, incx, शून्य, 0, शून्य, 1);

#अगर_घोषित SMP
  } अन्यथा {

#अगर_घोषित DOUBLE
    mode  =  BLAS_DOUBLE | BLAS_REAL;
#अन्यथा
    mode  =  BLAS_SINGLE | BLAS_REAL;
#पूर्ण_अगर

    blas_level1_thread(mode, n, 0, 0,
#अगर_अघोषित CBLAS
		       ALPHA,
#अन्यथा
		       &alpha,
#पूर्ण_अगर
		       x, incx, शून्य, 0, शून्य, 1,  (पूर्णांक (*)(व्योम))SCAL_K, nthreads);

  }
#पूर्ण_अगर

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  वापस;

}
