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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XSPR  "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSPR  "
#มิฉะนั้น
#define ERROR_NAME "CSPR  "
#endif

static จำนวนเต็ม (*spr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xspr_U, xspr_L,
#elif defined(DOUBLE)
  zspr_U, zspr_L,
#มิฉะนั้น
  cspr_U, cspr_L,
#endif
};

#ifdef SMP
static จำนวนเต็ม (*spr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, FLOAT *, จำนวนเต็ม) = {
#ifdef XDOUBLE
  xspr_thread_U, xspr_thread_L,
#elif defined(DOUBLE)
  zspr_thread_U, zspr_thread_L,
#มิฉะนั้น
  cspr_thread_U, cspr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;

  blasint info;
  จำนวนเต็ม uplo;
  FLOAT *buffer;
#ifdef SMP
  จำนวนเต็ม nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  ถ้า (uplo_arg  == 'U') uplo  = 0;
  ถ้า (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ถ้า (incx == 0)          info =  5;
  ถ้า (n < 0)              info =  2;
  ถ้า (uplo  < 0)          info =  1;

  ถ้า (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

  ถ้า (n == 0) คืน;

  ถ้า ((alpha_r == ZERO) && (alpha_i == ZERO)) คืน;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ถ้า (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ถ้า (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    คืน;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ถ้า (n < 50) nthreads = 1;
  ถ้า (nthreads > 2 && n < 500) nthreads = 2;

  ถ้า (nthreads == 1) {
#endif

    (spr[uplo])(n, alpha_r, alpha_i, x, incx, a, buffer);

#ifdef SMP
  } มิฉะนั้น {

    (spr_thread[uplo])(n, ALPHA, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  คืน;
}
