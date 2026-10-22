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
#define ERROR_NAME "QSPMV "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DSPMV "
#अन्यथा
#define ERROR_NAME "SSPMV "
#endif

स्थिर  पूर्णांक (*spmv[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  qspmv_U, qspmv_L,
#किंवाजर defined(DOUBLE)
  dspmv_U, dspmv_L,
#अन्यथा
  sspmv_U, sspmv_L,
#endif
};

#ifdef SMPTEST
स्थिर  पूर्णांक (*spmv_thread[])(BLASLONG, FLOAT, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#ifdef XDOUBLE
  qspmv_thread_U, qspmv_thread_L,
#किंवाजर defined(DOUBLE)
  dspmv_thread_U, dspmv_thread_L,
#अन्यथा
  sspmv_thread_U, sspmv_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *UPLO, blasint *N, FLOAT  *ALPHA, FLOAT *a,
            FLOAT  *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  FLOAT beta	= *BETA;
  blasint incy	= *INCY;

  blasint info;
  पूर्णांक uplo;
  FLOAT *buffer;
#ifdef SMPTEST
  पूर्णांक nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  जर (uplo_arg  == 'U') uplo  = 0;
  जर (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  जर (incy == 0)          info =  9;
  जर (incx == 0)          info =  6;
  जर (n < 0)              info =  2;
  जर (uplo  < 0)          info =  1;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *a,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  पूर्णांक uplo;
  blasint info;
#ifdef SMPTEST
  पूर्णांक nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  जर (order == CblasColMajor) {
    जर (Uplo == CblasUpper)         uplo  = 0;
    जर (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    जर (incy == 0)          info =  9;
    जर (incx == 0)          info =  6;
    जर (n < 0)              info =  2;
    जर (uplo  < 0)          info =  1;
  }

  जर (order == CblasRowMajor) {
    जर (Uplo == CblasUpper)         uplo  = 1;
    जर (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    जर (incy == 0)          info =  9;
    जर (incx == 0)          info =  6;
    जर (n < 0)              info =  2;
    जर (uplo  < 0)          info =  1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#endif

  जर (n == 0) परत;

  जर (beta != ONE) SCAL_K(n, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  जर (alpha == ZERO) परत;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जर (incx < 0 ) x -= (n - 1) * incx;
  जर (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#ifdef SMPTEST
  nthreads = अंक_cpu_avail(2);

  जर (nthreads == 1) {
#endif

  (spmv[uplo])(n, alpha, a, x, incx, y, incy, buffer);

#ifdef SMPTEST
  } अन्यथा {

    (spmv_thread[uplo])(n, alpha, a, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  परत;
}
