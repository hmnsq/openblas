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

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTPMV "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DTPMV "
#وگرنه
#define ERROR_NAME "STPMV "
#endif

ایستا صحیح (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, تهی *) = {
#ifdef XDOUBLE
  qtpmv_NUU, qtpmv_NUN, qtpmv_NLU, qtpmv_NLN,
  qtpmv_TUU, qtpmv_TUN, qtpmv_TLU, qtpmv_TLN,
#وگرنهاگر defined(DOUBLE)
  dtpmv_NUU, dtpmv_NUN, dtpmv_NLU, dtpmv_NLN,
  dtpmv_TUU, dtpmv_TUN, dtpmv_TLU, dtpmv_TLN,
#وگرنه
  stpmv_NUU, stpmv_NUN, stpmv_NLU, stpmv_NLN,
  stpmv_TUU, stpmv_TUN, stpmv_TLU, stpmv_TLN,
#endif
};

#ifdef SMP
ایستا صحیح (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#ifdef XDOUBLE
  qtpmv_thread_NUU, qtpmv_thread_NUN, qtpmv_thread_NLU, qtpmv_thread_NLN,
  qtpmv_thread_TUU, qtpmv_thread_TUN, qtpmv_thread_TLU, qtpmv_thread_TLN,
#وگرنهاگر defined(DOUBLE)
  dtpmv_thread_NUU, dtpmv_thread_NUN, dtpmv_thread_NLU, dtpmv_thread_NLN,
  dtpmv_thread_TUU, dtpmv_thread_TUN, dtpmv_thread_TLU, dtpmv_thread_TLN,
#وگرنه
  stpmv_thread_NUU, stpmv_thread_NUN, stpmv_thread_NLU, stpmv_thread_NLN,
  stpmv_thread_TUU, stpmv_thread_TUN, stpmv_thread_TLU, stpmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(نویسه *UPLO, نویسه *TRANS, نویسه *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  نویسه uplo_arg  = *UPLO;
  نویسه trans_arg = *TRANS;
  نویسه diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  صحیح uplo;
  صحیح unit;
  صحیح trans;
  FLOAT *buffer;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'R') trans = 0;
  اگر (trans_arg == 'C') trans = 1;

  اگر (diag_arg  == 'U') unit  = 0;
  اگر (diag_arg  == 'N') unit  = 1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  اگر (incx == 0)          info =  7;
  اگر (n < 0)              info =  4;
  اگر (unit  < 0)          info =  3;
  اگر (trans < 0)          info =  2;
  اگر (uplo  < 0)          info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#وگرنه
#ifndef COMPLEX
تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_UPLO Uplo,
	   شمارش CBLAS_TRANSPOSE TransA, شمارش CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {
#وگرنه
تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_UPLO Uplo,
	   شمارش CBLAS_TRANSPOSE TransA, شمارش CBLAS_DIAG Diag,
	   blasint n, تهی  *va, تهی  *vx, blasint incx) {
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
#endif

  صحیح trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  صحیح nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  اگر (order == CblasColMajor) {
    اگر (Uplo == CblasUpper)         uplo  = 0;
    اگر (Uplo == CblasLower)         uplo  = 1;

    اگر (TransA == CblasNoTrans)     trans = 0;
    اگر (TransA == CblasTrans)       trans = 1;
    اگر (TransA == CblasConjNoTrans) trans = 0;
    اگر (TransA == CblasConjTrans)   trans = 1;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  7;
    اگر (n < 0)              info =  4;
    اگر (unit  < 0)          info =  3;
    اگر (trans < 0)          info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (order == CblasRowMajor) {
    اگر (Uplo == CblasUpper)         uplo  = 1;
    اگر (Uplo == CblasLower)         uplo  = 0;

    اگر (TransA == CblasNoTrans)     trans = 1;
    اگر (TransA == CblasTrans)       trans = 0;
    اگر (TransA == CblasConjNoTrans) trans = 1;
    اگر (TransA == CblasConjTrans)   trans = 0;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  7;
    اگر (n < 0)              info =  4;
    اگر (unit  < 0)          info =  3;
    اگر (trans < 0)          info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#endif

  اگر (n == 0) بازگشت;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  اگر (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  اگر (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
    بازگشت;
  }

#ifdef SMP
  nthreads = رقم_cpu_avail(2);

  اگر (n < 50) nthreads = 1;
  اگر (nthreads > 2 && n < 500) nthreads = 2;

  اگر (nthreads == 1) {
#endif

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP
  } وگرنه {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  بازگشت;
}
