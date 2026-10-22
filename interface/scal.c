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
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#وإلا

فراغ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#جے_انت

#جے_تعریف_شدہ SMP
  صحيح منوال, nthreads;
#جے_انت

#جے_بنا_تعریف CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جے_انت

  إذا (incx <= 0 || n <= 0) أرجع;

#جے_تعریف_شدہ BGEMM
  عائم alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#وإلا
  FLOAT alpha_float = alpha;
#جے_انت

  إذا (alpha_float == ONE) أرجع;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();


#جے_تعریف_شدہ SMP
  إذا (n <= 1048576 )
	nthreads = 1;
  وإلا
	nthreads = عدد_cpu_avail(1);


  إذا (nthreads == 1) {
#جے_انت

  SCAL_K(n, 0, 0, alpha, x, incx, عدم, 0, عدم, 1);

#جے_تعریف_شدہ SMP
  } وإلا {

#جے_تعریف_شدہ DOUBLE
    منوال  =  BLAS_DOUBLE | BLAS_REAL;
#وإلا
    منوال  =  BLAS_SINGLE | BLAS_REAL;
#جے_انت

    blas_level1_thread(منوال, n, 0, 0,
#جے_بنا_تعریف CBLAS
		       ALPHA,
#وإلا
		       &alpha,
#جے_انت
		       x, incx, عدم, 0, عدم, 1,  (صحيح (*)(فراغ))SCAL_K, nthreads);

  }
#جے_انت

  FUNCTION_PROFILE_سرا(1, n, n);

  IDEBUG_سرا;

  أرجع;

}
