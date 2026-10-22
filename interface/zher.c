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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XHER  "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZHER  "
#अन्यथा
#परिभाषय ERROR_NAME "CHER  "
#यदि_अन्तः

स्थिर पूर्णाङ्क (*her[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#यदि_परिभाषितम् XDOUBLE
  xher_U, xher_L, xher_V, xher_M,
#नोचेद्यदि defined(DOUBLE)
  zher_U, zher_L, zher_V, zher_M,
#अन्यथा
  cher_U, cher_L, cher_V, cher_M,
#यदि_अन्तः
};

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*her_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#यदि_परिभाषितम् XDOUBLE
  xher_thread_U, xher_thread_L, xher_thread_V, xher_thread_M,
#नोचेद्यदि defined(DOUBLE)
  zher_thread_U, zher_thread_L, zher_thread_V, zher_thread_M,
#अन्यथा
  cher_thread_U, cher_thread_L, cher_thread_V, cher_thread_M,
#यदि_अन्तः
};
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	  FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *बफर;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  यदि (lda  < MAX(1, n))   info =  7;
  यदि (incx == 0)          info =  5;
  यदि (n < 0)              info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo, blasint n, FLOAT alpha, शून्य *vx, blasint incx, शून्य *va, blasint lda) {

  FLOAT* a = (FLOAT*) va;
  FLOAT* x = (FLOAT*) vx;

  FLOAT *बफर;
  पूर्णाङ्क uplo;
  blasint info;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 1;

    info = -1;

    यदि (lda  < MAX(1, n))   info =  7;
    यदि (incx == 0)          info =  5;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;

  }

  यदि (क्रमः == CblasRowMajor) {

    यदि (Uplo == Cblas_ऊर्ध्वस्थम्) uplo  = 3;
    यदि (Uplo == Cblas_अधस्तनम्) uplo  = 2;

    info = -1;

    यदि (lda  < MAX(1, n))   info =  7;
    यदि (incx == 0)          info =  5;
    यदि (n < 0)              info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  यदि (n == 0) प्रत्यागम;

  यदि (alpha == ZERO) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx * 2;

  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_परिभाषितम् SMP
  यदि (n < 100)
    nthreads = 1;
  अन्यथा
    nthreads = अङ्कः_cpu_avail(2);

  यदि (nthreads == 1) {
#यदि_अन्तः

  (her[uplo])(n, alpha, x, incx, a, lda, बफर);

#यदि_परिभाषितम् SMP
  } अन्यथा {

    (her_thread[uplo])(n, alpha, x, incx, a, lda, बफर, nthreads);

  }
#यदि_अन्तः

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(4, n * n / 2 + n, n * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
