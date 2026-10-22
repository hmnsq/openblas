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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

فراغ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#وإلا

فراغ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#د_که_پای

#که_تعریف_شوی SMP
  صحيح ډېرپېښ, nthreads;
#د_که_پای

#که_نه_تعریف_شوی CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#د_که_پای

  إذا (incx <= 0 || n <= 0) أرجع;

#که_تعریف_شوی BGEMM
  عائم alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#وإلا
  FLOAT alpha_float = alpha;
#د_که_پای

  إذا (alpha_float == ONE) أرجع;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();


#که_تعریف_شوی SMP
  إذا (n <= 1048576 )
	nthreads = 1;
  وإلا
	nthreads = عدد_cpu_avail(1);


  إذا (nthreads == 1) {
#د_که_پای

  SCAL_K(n, 0, 0, alpha, x, incx, عدم, 0, عدم, 1);

#که_تعریف_شوی SMP
  } وإلا {

#که_تعریف_شوی DOUBLE
    ډېرپېښ  =  BLAS_DOUBLE | BLAS_REAL;
#وإلا
    ډېرپېښ  =  BLAS_SINGLE | BLAS_REAL;
#د_که_پای

    blas_level1_thread(ډېرپېښ, n, 0, 0,
#که_نه_تعریف_شوی CBLAS
		       ALPHA,
#وإلا
		       &alpha,
#د_که_پای
		       x, incx, عدم, 0, عدم, 1,  (صحيح (*)(فراغ))SCAL_K, nthreads);

  }
#د_که_پای

  FUNCTION_PROFILE_پای_ټکی(1, n, n);

  IDEBUG_پای_ټکی;

  أرجع;

}
