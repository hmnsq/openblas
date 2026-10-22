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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <ctype.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QTPMV "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DTPMV "
#இல்லையெனில்
#define ERROR_NAME "STPMV "
#endif

நிலையான முழு (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, வெறுமை *) = {
#ifdef XDOUBLE
  qtpmv_NUU, qtpmv_NUN, qtpmv_NLU, qtpmv_NLN,
  qtpmv_TUU, qtpmv_TUN, qtpmv_TLU, qtpmv_TLN,
#அல்லதெனில் defined(DOUBLE)
  dtpmv_NUU, dtpmv_NUN, dtpmv_NLU, dtpmv_NLN,
  dtpmv_TUU, dtpmv_TUN, dtpmv_TLU, dtpmv_TLN,
#இல்லையெனில்
  stpmv_NUU, stpmv_NUN, stpmv_NLU, stpmv_NLN,
  stpmv_TUU, stpmv_TUN, stpmv_TLU, stpmv_TLN,
#endif
};

#ifdef SMP
நிலையான முழு (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, முழு) = {
#ifdef XDOUBLE
  qtpmv_thread_NUU, qtpmv_thread_NUN, qtpmv_thread_NLU, qtpmv_thread_NLN,
  qtpmv_thread_TUU, qtpmv_thread_TUN, qtpmv_thread_TLU, qtpmv_thread_TLN,
#அல்லதெனில் defined(DOUBLE)
  dtpmv_thread_NUU, dtpmv_thread_NUN, dtpmv_thread_NLU, dtpmv_thread_NLN,
  dtpmv_thread_TUU, dtpmv_thread_TUN, dtpmv_thread_TLU, dtpmv_thread_TLN,
#இல்லையெனில்
  stpmv_thread_NUU, stpmv_thread_NUN, stpmv_thread_NLU, stpmv_thread_NLN,
  stpmv_thread_TUU, stpmv_thread_TUN, stpmv_thread_TLU, stpmv_thread_TLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, எழுத்து *TRANS, எழுத்து *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  எழுத்து uplo_arg  = *UPLO;
  எழுத்து trans_arg = *TRANS;
  எழுத்து diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  முழு uplo;
  முழு unit;
  முழு trans;
  FLOAT *buffer;
#ifdef SMP
  முழு nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
  எனில் (trans_arg == 'R') trans = 0;
  எனில் (trans_arg == 'C') trans = 1;

  எனில் (diag_arg  == 'U') unit  = 0;
  எனில் (diag_arg  == 'N') unit  = 1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  எனில் (incx == 0)          info =  7;
  எனில் (n < 0)              info =  4;
  எனில் (unit  < 0)          info =  3;
  எனில் (trans < 0)          info =  2;
  எனில் (uplo  < 0)          info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்
#ifndef COMPLEX
வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA, எண்ணிக்கை CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {
#இல்லையெனில்
வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo,
	   எண்ணிக்கை CBLAS_TRANSPOSE TransA, எண்ணிக்கை CBLAS_DIAG Diag,
	   blasint n, வெறுமை  *va, வெறுமை  *vx, blasint incx) {
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
#endif

  முழு trans, uplo, unit;
  blasint info;
  FLOAT *buffer;
#ifdef SMP
  முழு nthreads;
#endif

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 0;
    எனில் (Uplo == CblasLower)         uplo  = 1;

    எனில் (TransA == CblasNoTrans)     trans = 0;
    எனில் (TransA == CblasTrans)       trans = 1;
    எனில் (TransA == CblasConjNoTrans) trans = 0;
    எனில் (TransA == CblasConjTrans)   trans = 1;

    எனில் (Diag == CblasUnit)          unit  = 0;
    எனில் (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    எனில் (incx == 0)          info =  7;
    எனில் (n < 0)              info =  4;
    எனில் (unit  < 0)          info =  3;
    எனில் (trans < 0)          info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (Uplo == CblasUpper)         uplo  = 1;
    எனில் (Uplo == CblasLower)         uplo  = 0;

    எனில் (TransA == CblasNoTrans)     trans = 1;
    எனில் (TransA == CblasTrans)       trans = 0;
    எனில் (TransA == CblasConjNoTrans) trans = 1;
    எனில் (TransA == CblasConjTrans)   trans = 0;

    எனில் (Diag == CblasUnit)          unit  = 0;
    எனில் (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    எனில் (incx == 0)          info =  7;
    எனில் (n < 0)              info =  4;
    எனில் (unit  < 0)          info =  3;
    எனில் (trans < 0)          info =  2;
    எனில் (uplo  < 0)          info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#endif

  எனில் (n == 0) திரும்பு;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  எனில் (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#ifdef SMP
  nthreads = இலக்கம்_cpu_avail(2);

  எனில் (n < 50) nthreads = 1;
  எனில் (nthreads > 2 && n < 500) nthreads = 2;

  எனில் (nthreads == 1) {
#endif

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP
  } இல்லையெனில் {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  திரும்பு;
}
