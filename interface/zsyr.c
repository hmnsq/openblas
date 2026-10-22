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
#define ERROR_NAME "XSYR  "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZSYR  "
#अन्यथा
#define ERROR_NAME "CSYR  "
#endif

स्थिर पूर्णांक (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr_U, xsyr_L,
#किंवाजर defined(DOUBLE)
  zsyr_U, zsyr_L,
#अन्यथा
  csyr_U, csyr_L,
#endif
};

#ifdef SMP
स्थिर पूर्णांक (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#ifdef XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#किंवाजर defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#अन्यथा
  csyr_thread_U, csyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(अक्षर *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  अक्षर uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  पूर्णांक uplo;
  FLOAT *buffer;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  जर (uplo_arg  == 'U') uplo  = 0;
  जर (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  जर (lda  < MAX(1, n))   info =  7;
  जर (incx == 0)          info =  5;
  जर (n < 0)              info =  2;
  जर (uplo  < 0)          info =  1;

  जर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }


#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo, पूर्णांक n, शून्य* valpha, FLOAT *x, पूर्णांक incx, FLOAT *a, पूर्णांक lda) {

  FLOAT *buffer;
  पूर्णांक uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
#ifdef SMP
  पूर्णांक nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  जर (order == CblasColMajor) {

    जर (Uplo == CblasUpper) uplo  = 0;
    जर (Uplo == CblasLower) uplo  = 1;

    info = -1;

    जर (lda  < MAX(1, n))   info =  7;
    जर (incx == 0)          info =  5;
    जर (n < 0)              info =  2;
    जर (uplo  < 0)          info =  1;

  }

  जर (order == CblasRowMajor) {

    जर (Uplo == CblasUpper) uplo  = 1;
    जर (Uplo == CblasLower) uplo  = 0;

    info = -1;

    जर (lda  < MAX(1, n))   info =  7;
    जर (incx == 0)          info =  5;
    जर (n < 0)              info =  2;
    जर (uplo  < 0)          info =  1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#endif

  जर (n == 0) परत;

  जर ((alpha_r == ZERO) && (alpha_i == ZERO)) परत;

  जर (incx == 1 && n < 50) {
    blasint i;
    जर (!uplo) {
      साठी (i = 0; i < n; i++){
        जर ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, रिक्त, 0);
        }
        a += lda * 2;
      }
    } अन्यथा {
      साठी (i = 0; i < n; i++){
        जर ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, रिक्त, 0);
        }
        a += 2 + lda * 2;
      }
    }
    परत;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जर (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#ifdef SMP
  nthreads = अंक_cpu_avail(2);

  जर (n < 50) nthreads = 1;
  जर (nthreads > 2 && n < 500) nthreads = 2;

  जर (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, buffer);

#ifdef SMP
  } अन्यथा {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  परत;
}
