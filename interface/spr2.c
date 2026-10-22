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
#define ERROR_NAME "QSPR2 "
#अथवायदि defined(DOUBLE)
#define ERROR_NAME "DSPR2 "
#अन्यथा
#define ERROR_NAME "SSPR2 "
#endif

स्थिर पूर्णाङ्क (*spr2[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr2_U, qspr2_L,
#अथवायदि defined(DOUBLE)
  dspr2_U, dspr2_L,
#अन्यथा
  sspr2_U, sspr2_L,
#endif
};

#ifdef SMP
स्थिर पूर्णाङ्क (*spr2_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, FLOAT *, पूर्णाङ्क) = {
#ifdef XDOUBLE
  qspr2_thread_U, qspr2_thread_L,
#अथवायदि defined(DOUBLE)
  dspr2_thread_U, dspr2_thread_L,
#अन्यथा
  sspr2_thread_U, sspr2_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *a){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;
  blasint incy	= *INCY;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *buffer;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  यदि (incy == 0)          info =  7;
  यदि (incx == 0)          info =  5;
  यदि (n < 0)              info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *y, blasint incy,
	   FLOAT  *a) {

  FLOAT *buffer;
  पूर्णाङ्क uplo;
  blasint info;
#ifdef SMP
  पूर्णाङ्क nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  यदि (order == CblasColMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 0;
    यदि (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    यदि (incy == 0)          info =  7;
    यदि (incx == 0)          info =  5;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (order == CblasRowMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 1;
    यदि (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    यदि (incy == 0)          info =  7;
    यदि (incx == 0)          info =  5;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#endif

  यदि (n == 0) फर्काऊ;

  यदि (alpha == ZERO) फर्काऊ;

  यदि (incx == 1 && incy == 1 && n < 50) {
    blasint i;
    यदि (!uplo) {
      लागि (i = 0; i < n; i++){
        AXPYU_K(i + 1, 0, 0, alpha * x[i], y,     1, a, 1, रिक्त, 0);
        AXPYU_K(i + 1, 0, 0, alpha * y[i], x,     1, a, 1, रिक्त, 0);
        a += i + 1;
      }
    } अन्यथा {
      लागि (i = 0; i < n; i++){
	AXPYU_K(n - i, 0, 0, alpha * x[i], y + i, 1, a, 1, रिक्त, 0);
        AXPYU_K(n - i, 0, 0, alpha * y[i], x + i, 1, a, 1, रिक्त, 0);
        a += n - i;
      }
    }
    फर्काऊ;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx;
  यदि (incy < 0 ) y -= (n - 1) * incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  यदि (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

#ifdef SMP
  nthreads = अङ्क_cpu_avail(2);

  यदि (n < 50) nthreads = 1;
  यदि (nthreads > 2 && n < 500) nthreads = 2;

  यदि (nthreads == 1) {
#endif

  (spr2[uplo])(n, alpha, x, incx, y, incy, a, buffer);

#ifdef SMP
  } अन्यथा {

    (spr2_thread[uplo])(n, alpha, x, incx, y, incy, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + 2 * n, 2 * n * n);

  IDEBUG_END;

  फर्काऊ;
}
