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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QTRSV "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DTRSV "
#وإلا
#وصف_ڪريو ERROR_NAME "STRSV "
#جيڪڏهن_پورو

ساكن صحيح (*trsv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, فراغ *) = {
#جيڪڏهن_وصف_ٿيل XDOUBLE
  qtrsv_NUU, qtrsv_NUN, qtrsv_NLU, qtrsv_NLN,
  qtrsv_TUU, qtrsv_TUN, qtrsv_TLU, qtrsv_TLN,
#وإلاإذا defined(DOUBLE)
  dtrsv_NUU, dtrsv_NUN, dtrsv_NLU, dtrsv_NLN,
  dtrsv_TUU, dtrsv_TUN, dtrsv_TLU, dtrsv_TLN,
#وإلا
  strsv_NUU, strsv_NUN, strsv_NLU, strsv_NLN,
  strsv_TUU, strsv_TUN, strsv_TLU, strsv_TLN,
#جيڪڏهن_پورو
};

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  صحيح uplo;
  صحيح unit;
  صحيح trans;
  FLOAT *بفر;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'R') trans = 0;
  إذا (trans_arg == 'C') trans = 1;

  إذا (diag_arg  == 'U') unit  = 0;
  إذا (diag_arg  == 'N') unit  = 1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  إذا (incx == 0)          info =  8;
  إذا (lda  < MAX(1, n))   info =  6;
  إذا (n < 0)              info =  4;
  إذا (unit  < 0)          info =  3;
  إذا (trans < 0)          info =  2;
  إذا (uplo  < 0)          info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_UPLO Uplo,
	   تعداد CBLAS_TRANSPOSE TransA, تعداد CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  صحيح trans, uplo, unit;
  blasint info;
  FLOAT *بفر;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (Uplo == Cblas_مٿانهون)         uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون)         uplo  = 1;

    إذا (TransA == CblasNoTrans)     trans = 0;
    إذا (TransA == CblasTrans)       trans = 1;
    إذا (TransA == CblasConjNoTrans) trans = 0;
    إذا (TransA == CblasConjTrans)   trans = 1;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  8;
    إذا (lda  < MAX(1, n))   info =  6;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (Uplo == Cblas_مٿانهون)         uplo  = 1;
    إذا (Uplo == Cblas_هيٺاهون)         uplo  = 0;

    إذا (TransA == CblasNoTrans)     trans = 1;
    إذا (TransA == CblasTrans)       trans = 0;
    إذا (TransA == CblasConjNoTrans) trans = 1;
    إذا (TransA == CblasConjTrans)   trans = 0;

    إذا (Diag == CblasUnit)          unit  = 0;
    إذا (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    إذا (incx == 0)          info =  8;
    إذا (lda  < MAX(1, n))   info =  6;
    إذا (n < 0)              info =  4;
    إذا (unit  < 0)          info =  3;
    إذا (trans < 0)          info =  2;
    إذا (uplo  < 0)          info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#جيڪڏهن_پورو

  إذا (n == 0) أرجع;

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  إذا (incx < 0 ) x -= (n - 1) * incx;

  بفر = (FLOAT *)blas_يادگيري_alloc(1);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  (trsv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, بفر);

  blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(1, n * n / 2 + n,  n * n);

  IDEBUG_پڇاڙي;

  أرجع;
}
