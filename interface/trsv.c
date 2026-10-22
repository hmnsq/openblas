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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <ctype.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTRSV "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DTRSV "
#નહીંતર
#define ERROR_NAME "STRSV "
#endif

સ્થિર પૂર્ણાંક (*trsv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, શૂન્ય *) = {
#ifdef XDOUBLE
  qtrsv_NUU, qtrsv_NUN, qtrsv_NLU, qtrsv_NLN,
  qtrsv_TUU, qtrsv_TUN, qtrsv_TLU, qtrsv_TLN,
#અથવાજો defined(DOUBLE)
  dtrsv_NUU, dtrsv_NUN, dtrsv_NLU, dtrsv_NLN,
  dtrsv_TUU, dtrsv_TUN, dtrsv_TLU, dtrsv_TLN,
#નહીંતર
  strsv_NUU, strsv_NUN, strsv_NLU, strsv_NLN,
  strsv_TUU, strsv_TUN, strsv_TLU, strsv_TLN,
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, અક્ષર *TRANS, અક્ષર *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  અક્ષર uplo_arg  = *UPLO;
  અક્ષર trans_arg = *TRANS;
  અક્ષર diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  પૂર્ણાંક uplo;
  પૂર્ણાંક unit;
  પૂર્ણાંક trans;
  FLOAT *buffer;

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'T') trans = 1;
  જો (trans_arg == 'R') trans = 0;
  જો (trans_arg == 'C') trans = 1;

  જો (diag_arg  == 'U') unit  = 0;
  જો (diag_arg  == 'N') unit  = 1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  જો (incx == 0)          info =  8;
  જો (lda  < MAX(1, n))   info =  6;
  જો (n < 0)              info =  4;
  જો (unit  < 0)          info =  3;
  જો (trans < 0)          info =  2;
  જો (uplo  < 0)          info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_UPLO Uplo,
	   ગણના CBLAS_TRANSPOSE TransA, ગણના CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, blasint lda, FLOAT  *x, blasint incx) {

  પૂર્ણાંક trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (Uplo == CblasUpper)         uplo  = 0;
    જો (Uplo == CblasLower)         uplo  = 1;

    જો (TransA == CblasNoTrans)     trans = 0;
    જો (TransA == CblasTrans)       trans = 1;
    જો (TransA == CblasConjNoTrans) trans = 0;
    જો (TransA == CblasConjTrans)   trans = 1;

    જો (Diag == CblasUnit)          unit  = 0;
    જો (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    જો (incx == 0)          info =  8;
    જો (lda  < MAX(1, n))   info =  6;
    જો (n < 0)              info =  4;
    જો (unit  < 0)          info =  3;
    જો (trans < 0)          info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Uplo == CblasUpper)         uplo  = 1;
    જો (Uplo == CblasLower)         uplo  = 0;

    જો (TransA == CblasNoTrans)     trans = 1;
    જો (TransA == CblasTrans)       trans = 0;
    જો (TransA == CblasConjNoTrans) trans = 1;
    જો (TransA == CblasConjTrans)   trans = 0;

    જો (Diag == CblasUnit)          unit  = 0;
    જો (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    જો (incx == 0)          info =  8;
    જો (lda  < MAX(1, n))   info =  6;
    જો (n < 0)              info =  4;
    જો (unit  < 0)          info =  3;
    જો (trans < 0)          info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  જો (n == 0) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  (trsv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, buffer);

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n,  n * n);

  IDEBUG_END;

  પાછા;
}
