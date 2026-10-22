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
#ضمن <ctype.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XSPR  "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSPR  "
#وإلا
#define ERROR_NAME "CSPR  "
#endif

ساكن صحيح (*spr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xspr_U, xspr_L,
#وإلاإذا defined(DOUBLE)
  zspr_U, zspr_L,
#وإلا
  cspr_U, cspr_L,
#endif
};

#ifdef SMP
ساكن صحيح (*spr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, FLOAT *, صحيح) = {
#ifdef XDOUBLE
  xspr_thread_U, xspr_thread_L,
#وإلاإذا defined(DOUBLE)
  zspr_thread_U, zspr_thread_L,
#وإلا
  cspr_thread_U, cspr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT
فراغ NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;

  blasint info;
  صحيح uplo;
  FLOAT *buffer;
#ifdef SMP
  صحيح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  𐤀𐤌 (uplo_arg  == 'U') uplo  = 0;
  𐤀𐤌 (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  𐤀𐤌 (incx == 0)          info =  5;
  𐤀𐤌 (n < 0)              info =  2;
  𐤀𐤌 (uplo  < 0)          info =  1;

  𐤀𐤌 (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  𐤀𐤌 (n == 0) أرجع;

  𐤀𐤌 ((alpha_r == ZERO) && (alpha_i == ZERO)) أرجع;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  𐤀𐤌 (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;

  buffer = (FLOAT *)blas_memory_alloc(1);
  𐤀𐤌 (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  𐤀𐤌 (n < 50) nthreads = 1;
  𐤀𐤌 (nthreads > 2 && n < 500) nthreads = 2;

  𐤀𐤌 (nthreads == 1) {
#endif

    (spr[uplo])(n, alpha_r, alpha_i, x, incx, a, buffer);

#ifdef SMP
  } وإلا {

    (spr_thread[uplo])(n, ALPHA, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  أرجع;
}
