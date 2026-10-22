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
#define ERROR_NAME "QTPSV "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DTPSV "
#अण्णहा
#define ERROR_NAME "STPSV "
#endif

स्थिर पूर्णाङ्क (*tpsv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  qtpsv_NUU, qtpsv_NUN, qtpsv_NLU, qtpsv_NLN,
  qtpsv_TUU, qtpsv_TUN, qtpsv_TLU, qtpsv_TLN,
#नोचेद्यदि defined(DOUBLE)
  dtpsv_NUU, dtpsv_NUN, dtpsv_NLU, dtpsv_NLN,
  dtpsv_TUU, dtpsv_TUN, dtpsv_TLU, dtpsv_TLN,
#अण्णहा
  stpsv_NUU, stpsv_NUN, stpsv_NLU, stpsv_NLN,
  stpsv_TUU, stpsv_TUN, stpsv_TLU, stpsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  जइ (trans_arg == 'N') trans = 0;
  जइ (trans_arg == 'T') trans = 1;
  जइ (trans_arg == 'R') trans = 0;
  जइ (trans_arg == 'C') trans = 1;

  जइ (diag_arg  == 'U') unit  = 0;
  जइ (diag_arg  == 'N') unit  = 1;

  जइ (uplo_arg  == 'U') uplo  = 0;
  जइ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  जइ (incx == 0)          info =  7;
  जइ (n < 0)              info =  4;
  जइ (unit  < 0)          info =  3;
  जइ (trans < 0)          info =  2;
  जइ (uplo  < 0)          info =  1;

  जइ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अण्णहा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {

  पूर्णाङ्क trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  जइ (order == CblasColMajor) {
    जइ (Uplo == CblasUpper)         uplo  = 0;
    जइ (Uplo == CblasLower)         uplo  = 1;

    जइ (TransA == CblasNoTrans)     trans = 0;
    जइ (TransA == CblasTrans)       trans = 1;
    जइ (TransA == CblasConjNoTrans) trans = 0;
    जइ (TransA == CblasConjTrans)   trans = 1;

    जइ (Diag == CblasUnit)          unit  = 0;
    जइ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    जइ (incx == 0)          info =  7;
    जइ (n < 0)              info =  4;
    जइ (unit  < 0)          info =  3;
    जइ (trans < 0)          info =  2;
    जइ (uplo  < 0)          info =  1;
  }

  जइ (order == CblasRowMajor) {
    जइ (Uplo == CblasUpper)         uplo  = 1;
    जइ (Uplo == CblasLower)         uplo  = 0;

    जइ (TransA == CblasNoTrans)     trans = 1;
    जइ (TransA == CblasTrans)       trans = 0;
    जइ (TransA == CblasConjNoTrans) trans = 1;
    जइ (TransA == CblasConjTrans)   trans = 0;

    जइ (Diag == CblasUnit)          unit  = 0;
    जइ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    जइ (incx == 0)          info =  7;
    जइ (n < 0)              info =  4;
    जइ (unit  < 0)          info =  3;
    जइ (trans < 0)          info =  2;
    जइ (uplo  < 0)          info =  1;
  }

  जइ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#endif

  जइ (n == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जइ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  जइ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  (tpsv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  प्रत्यागम;
}
