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

#شامل <stdio.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;

#اگر_غیر_تعریف_شدہ SSCAL
  FLOAT *alpha=ALPHA;
#ورنہ
  FLOAT alpha[2] = {ALPHA[0], ZERO};
#اگر_ختم

#ورنہ

#اگر_غیر_تعریف_شدہ SSCAL
خالی CNAME(blasint n, خالی *VALPHA, خالی *vx, blasint incx){

  FLOAT *x = (FLOAT*) vx;
  FLOAT *alpha=(FLOAT*)VALPHA;
#ورنہ
خالی CNAME(blasint n, FLOAT alpha_r, خالی *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;

  FLOAT alpha[2] = {alpha_r, ZERO};
#اگر_ختم
#اگر_ختم

#اگر_تعریف_شدہ SMP
  صحیح منوال;
  صحیح nthreads;
#اگر_ختم

#اگر_غیر_تعریف_شدہ CBLAS
  PRINT_DEBUG_NAME;
#ورنہ
  PRINT_DEBUG_CNAME;
#اگر_ختم

  اگر (incx <= 0 || n <= 0) واپس;

  اگر ((alpha[0] == ONE) && (alpha[1] == ZERO)) واپس;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

#اگر_تعریف_شدہ SMP
  اگر ( n <= 1048576 )
	nthreads = 1;
  ورنہ
	nthreads = عدد_cpu_avail(1);

  اگر (nthreads == 1) {
#اگر_ختم

  SCAL_K(n, 0, 0, alpha[0], alpha[1], x, incx, تہی, 0, تہی, 1);

#اگر_تعریف_شدہ SMP
  } ورنہ {
#اگر_تعریف_شدہ DOUBLE
    منوال  =  BLAS_DOUBLE | BLAS_COMPLEX;
#ورنہ
    منوال  =  BLAS_SINGLE | BLAS_COMPLEX;
#اگر_ختم

    blas_level1_thread(منوال, n, 0, 0,  alpha, x, incx, تہی, 0, تہی, 1, (صحیح (*)(خالی))SCAL_K, nthreads);

  }
#اگر_ختم

  FUNCTION_PROFILE_انتہا(4, n, n);

  IDEBUG_انتہا;

  واپس;

}
