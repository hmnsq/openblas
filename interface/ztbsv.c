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
#define ERROR_NAME "XTBSV "
#अथवायदि defined(DOUBLE)
#define ERROR_NAME "ZTBSV "
#अन्यथा
#define ERROR_NAME "CTBSV "
#endif

स्थिर पूर्णाङ्क (*tbsv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
  xtbsv_NUU, xtbsv_NUN, xtbsv_NLU, xtbsv_NLN,
  xtbsv_TUU, xtbsv_TUN, xtbsv_TLU, xtbsv_TLN,
  xtbsv_RUU, xtbsv_RUN, xtbsv_RLU, xtbsv_RLN,
  xtbsv_CUU, xtbsv_CUN, xtbsv_CLU, xtbsv_CLN,
#अथवायदि defined(DOUBLE)
  ztbsv_NUU, ztbsv_NUN, ztbsv_NLU, ztbsv_NLN,
  ztbsv_TUU, ztbsv_TUN, ztbsv_TLU, ztbsv_TLN,
  ztbsv_RUU, ztbsv_RUN, ztbsv_RLU, ztbsv_RLN,
  ztbsv_CUU, ztbsv_CUN, ztbsv_CLU, ztbsv_CLN,
#अन्यथा
  ctbsv_NUU, ctbsv_NUN, ctbsv_NLU, ctbsv_NLN,
  ctbsv_TUU, ctbsv_TUN, ctbsv_TLU, ctbsv_TLN,
  ctbsv_RUU, ctbsv_RUN, ctbsv_RLU, ctbsv_RLN,
  ctbsv_CUU, ctbsv_CUN, ctbsv_CLU, ctbsv_CLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

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
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'R') trans = 2;
  यदि (trans_arg == 'C') trans = 3;

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
    फर्काऊ;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, blasint k, शून्य  *va, blasint lda, शून्य  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  पूर्णाङ्क trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  यदि (order == CblasColMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 0;
    यदि (Uplo == CblasLower)         uplo  = 1;

    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 2;
    यदि (TransA == CblasConjTrans)   trans = 3;

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

  यदि (order == CblasRowMajor) {
    यदि (Uplo == CblasUpper)         uplo  = 1;
    यदि (Uplo == CblasLower)         uplo  = 0;

    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 3;
    यदि (TransA == CblasConjTrans)   trans = 2;

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
    फर्काऊ;
  }

#endif

  यदि (n == 0) फर्काऊ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  यदि (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    फर्काऊ;
  }

  (tbsv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * k / 2 + n, n * k);

  IDEBUG_END;

  फर्काऊ;
}
