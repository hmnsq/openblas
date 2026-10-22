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

#ضمن <stdio.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#وإلا

فراغ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#پایان_اگر

#اگر_تعریف_شده SMP
  صحيح نما, nthreads;
#پایان_اگر

#اگر_تعریف_نشده CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#پایان_اگر

  إذا (incx <= 0 || n <= 0) أرجع;

#اگر_تعریف_شده BGEMM
  عائم alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#وإلا
  FLOAT alpha_float = alpha;
#پایان_اگر

  إذا (alpha_float == ONE) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();


#اگر_تعریف_شده SMP
  إذا (n <= 1048576 )
	nthreads = 1;
  وإلا
	nthreads = رقم_cpu_avail(1);


  إذا (nthreads == 1) {
#پایان_اگر

  SCAL_K(n, 0, 0, alpha, x, incx, عدم, 0, عدم, 1);

#اگر_تعریف_شده SMP
  } وإلا {

#اگر_تعریف_شده DOUBLE
    نما  =  BLAS_DOUBLE | BLAS_REAL;
#وإلا
    نما  =  BLAS_SINGLE | BLAS_REAL;
#پایان_اگر

    blas_level1_thread(نما, n, 0, 0,
#اگر_تعریف_نشده CBLAS
		       ALPHA,
#وإلا
		       &alpha,
#پایان_اگر
		       x, incx, عدم, 0, عدم, 1,  (صحيح (*)(فراغ))SCAL_K, nthreads);

  }
#پایان_اگر

  FUNCTION_PROFILE_انتها(1, n, n);

  IDEBUG_انتها;

  أرجع;

}
