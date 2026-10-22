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
#परिभाषय ERROR_NAME "QTBSV "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DTBSV "
#अन्यथा
#परिभाषय ERROR_NAME "STBSV "
#यदि_अन्तः

स्थिर पूर्णाङ्क (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#यदि_परिभाषितम् XDOUBLE
  qtbsv_NUU, qtbsv_NUN, qtbsv_NLU, qtbsv_NLN,
  qtbsv_TUU, qtbsv_TUN, qtbsv_TLU, qtbsv_TLN,
#नोचेद्यदि defined(DOUBLE)
  dtbsv_NUU, dtbsv_NUN, dtbsv_NLU, dtbsv_NLN,
  dtbsv_TUU, dtbsv_TUN, dtbsv_TLU, dtbsv_TLN,
#अन्यथा
  stbsv_NUU, stbsv_NUN, stbsv_NLU, stbsv_NLN,
  stbsv_TUU, stbsv_TUN, stbsv_TLU, stbsv_TLN,
#यदि_अन्तः
};

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans;
  FLOAT *बफर;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'R') trans = 0;
  यदि (trans_arg == 'C') trans = 1;

  यदि (diag_arg  == 'U') unit  = 0;
  यदि (diag_arg  == 'N') unit  = 1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  यदि (incx == 0)          info =  9;
  यदि (lda < k + 1)        info =  7;
  यदि (k < 0)              info =  5;
  यदि (n < 0)              info =  4;
  यदि (unit  < 0)          info =  3;
  यदि (trans < 0)          info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, blasint k, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  पूर्णाङ्क trans, uplo, unit;
  blasint info;
  FLOAT *बफर;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)         uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्)         uplo  = 1;

    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 0;
    यदि (TransA == CblasConjTrans)   trans = 1;

    यदि (Diag == CblasUnit)          unit  = 0;
    यदि (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    यदि (incx == 0)          info =  9;
    यदि (lda < k + 1)        info =  7;
    यदि (k < 0)              info =  5;
    यदि (n < 0)              info =  4;
    यदि (unit  < 0)          info =  3;
    यदि (trans < 0)          info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)         uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्)         uplo  = 0;

    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 1;
    यदि (TransA == CblasConjTrans)   trans = 0;

    यदि (Diag == CblasUnit)          unit  = 0;
    यदि (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    यदि (incx == 0)          info =  9;
    यदि (lda < k + 1)        info =  7;
    यदि (k < 0)              info =  5;
    यदि (n < 0)              info =  4;
    यदि (unit  < 0)          info =  3;
    यदि (trans < 0)          info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  यदि (n == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx;

  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, बफर);

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(1, n * k / 2 + n, n * k);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
