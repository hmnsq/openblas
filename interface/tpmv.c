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

#include <stdio.h>
#include <ctype.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTPMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DTPMV "
#ບໍ່ດັ່ງນັ້ນ
#define ERROR_NAME "STPMV "
#endif

static ຈຳນວນເຕັມ (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
  qtpmv_NUU, qtpmv_NUN, qtpmv_NLU, qtpmv_NLN,
  qtpmv_TUU, qtpmv_TUN, qtpmv_TLU, qtpmv_TLN,
#elif defined(DOUBLE)
  dtpmv_NUU, dtpmv_NUN, dtpmv_NLU, dtpmv_NLN,
  dtpmv_TUU, dtpmv_TUN, dtpmv_TLU, dtpmv_TLN,
#ບໍ່ດັ່ງນັ້ນ
  stpmv_NUU, stpmv_NUN, stpmv_NLU, stpmv_NLN,
  stpmv_TUU, stpmv_TUN, stpmv_TLU, stpmv_TLN,
#endif
};

#ifdef SMP
static ຈຳນວນເຕັມ (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, ຈຳນວນເຕັມ) = {
#ifdef XDOUBLE
  qtpmv_thread_NUU, qtpmv_thread_NUN, qtpmv_thread_NLU, qtpmv_thread_NLN,
  qtpmv_thread_TUU, qtpmv_thread_TUN, qtpmv_thread_TLU, qtpmv_thread_TLN,
#elif defined(DOUBLE)
  dtpmv_thread_NUU, dtpmv_thread_NUN, dtpmv_thread_NLU, dtpmv_thread_NLN,
  dtpmv_thread_TUU, dtpmv_thread_TUN, dtpmv_thread_TLU, dtpmv_thread_TLN,
#ບໍ່ດັ່ງນັ້ນ
  stpmv_thread_NUU, stpmv_thread_NUN, stpmv_thread_NLU, stpmv_thread_NLN,
  stpmv_thread_TUU, stpmv_thread_TUN, stpmv_thread_TLU, stpmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, var_na *TRANS, var_na *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  var_na uplo_arg  = *UPLO;
  var_na trans_arg = *TRANS;
  var_na diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  ຈຳນວນເຕັມ uplo;
  ຈຳນວນເຕັມ unit;
  ຈຳນວນເຕັມ trans;
  FLOAT *buffer;
#ifdef SMP
  ຈຳນວນເຕັມ nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  ຖ້າ (trans_arg == 'N') trans = 0;
  ຖ້າ (trans_arg == 'T') trans = 1;
  ຖ້າ (trans_arg == 'R') trans = 0;
  ຖ້າ (trans_arg == 'C') trans = 1;

  ຖ້າ (diag_arg  == 'U') unit  = 0;
  ຖ້າ (diag_arg  == 'N') unit  = 1;

  ຖ້າ (uplo_arg  == 'U') uplo  = 0;
  ຖ້າ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  ຖ້າ (incx == 0)          info =  7;
  ຖ້າ (n < 0)              info =  4;
  ຖ້າ (unit  < 0)          info =  3;
  ຖ້າ (trans < 0)          info =  2;
  ຖ້າ (uplo  < 0)          info =  1;

  ຖ້າ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#ບໍ່ດັ່ງນັ້ນ
#ifndef COMPLEX
void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {
#ບໍ່ດັ່ງນັ້ນ
void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo,
	   ga_nanaa CBLAS_TRANSPOSE TransA, ga_nanaa CBLAS_DIAG Diag,
	   blasint n, void  *va, void  *vx, blasint incx) {
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
#endif

  ຈຳນວນເຕັມ trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  ຈຳນວນເຕັມ nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  ຖ້າ (order == CblasColMajor) {
    ຖ້າ (Uplo == CblasUpper)         uplo  = 0;
    ຖ້າ (Uplo == CblasLower)         uplo  = 1;

    ຖ້າ (TransA == CblasNoTrans)     trans = 0;
    ຖ້າ (TransA == CblasTrans)       trans = 1;
    ຖ້າ (TransA == CblasConjNoTrans) trans = 0;
    ຖ້າ (TransA == CblasConjTrans)   trans = 1;

    ຖ້າ (Diag == CblasUnit)          unit  = 0;
    ຖ້າ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ຖ້າ (incx == 0)          info =  7;
    ຖ້າ (n < 0)              info =  4;
    ຖ້າ (unit  < 0)          info =  3;
    ຖ້າ (trans < 0)          info =  2;
    ຖ້າ (uplo  < 0)          info =  1;
  }

  ຖ້າ (order == CblasRowMajor) {
    ຖ້າ (Uplo == CblasUpper)         uplo  = 1;
    ຖ້າ (Uplo == CblasLower)         uplo  = 0;

    ຖ້າ (TransA == CblasNoTrans)     trans = 1;
    ຖ້າ (TransA == CblasTrans)       trans = 0;
    ຖ້າ (TransA == CblasConjNoTrans) trans = 1;
    ຖ້າ (TransA == CblasConjTrans)   trans = 0;

    ຖ້າ (Diag == CblasUnit)          unit  = 0;
    ຖ້າ (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    ຖ້າ (incx == 0)          info =  7;
    ຖ້າ (n < 0)              info =  4;
    ຖ້າ (unit  < 0)          info =  3;
    ຖ້າ (trans < 0)          info =  2;
    ຖ້າ (uplo  < 0)          info =  1;
  }

  ຖ້າ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#endif

  ຖ້າ (n == 0) ສົ່ງຄືນ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ຖ້າ (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  ຖ້າ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ສົ່ງຄືນ;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  ຖ້າ (n < 50) nthreads = 1;
  ຖ້າ (nthreads > 2 && n < 500) nthreads = 2;

  ຖ້າ (nthreads == 1) {
#endif

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP
  } ບໍ່ດັ່ງນັ້ນ {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  ສົ່ງຄືນ;
}
