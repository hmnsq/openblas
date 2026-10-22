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
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XHPR  "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZHPR  "
#अञ्ञथा
#define ERROR_NAME "CHPR  "
#endif

स्थिर पूर्णाङ्क (*hpr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  xhpr_U, xhpr_L, xhpr_V, xhpr_M,
#नोचेद्यदि defined(DOUBLE)
  zhpr_U, zhpr_L, zhpr_V, zhpr_M,
#अञ्ञथा
  chpr_U, chpr_L, chpr_V, chpr_M,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*hpr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  xhpr_thread_U, xhpr_thread_L, xhpr_thread_V, xhpr_thread_M,
#नोचेद्यदि defined(DOUBLE)
  zhpr_thread_U, zhpr_thread_L, zhpr_thread_V, zhpr_thread_M,
#अञ्ञथा
  chpr_thread_U, chpr_thread_L, chpr_thread_V, chpr_thread_M,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  सचे (uplo_arg  == 'U') uplo  = 0;
  सचे (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  सचे (incx == 0)          info =  5;
  सचे (n < 0)              info =  2;
  सचे (uplo  < 0)          info =  1;

  सचे (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   शून्य  *vx, blasint incx,
	   शून्य  *va) {

  FLOAT* x = (FLOAT*) vx;
  FLOAT* a = (FLOAT*) va;

  FLOAT *buffer;
  पूर्णाङ्क uplo;
  blasint info;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  सचे (order == CblasColMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 0;
    सचे (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    सचे (incx == 0)          info =  5;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (order == CblasRowMajor) {
    सचे (Uplo == CblasUpper)         uplo  = 3;
    सचे (Uplo == CblasLower)         uplo  = 2;

    info = -1;

    सचे (incx == 0)          info =  5;
    सचे (n < 0)              info =  2;
    सचे (uplo  < 0)          info =  1;
  }

  सचे (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#endif

  सचे (n == 0) पटिनिवत्त;

  सचे (alpha == ZERO) पटिनिवत्त;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  सचे (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  सचे (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    पटिनिवत्त;
  }

#ifdef SMP
  सचे (n < 100)
    nthreads = 1;
  अञ्ञथा
    nthreads = अङ्कः_cpu_avail(2);

  सचे (nthreads == 1) {
#endif

  (hpr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP

  } अञ्ञथा {

    (hpr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  पटिनिवत्त;
}
