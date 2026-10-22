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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "XTPMV "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZTPMV "
#অন্যথা
#define ERROR_NAME "CTPMV "
#endif

স্থির পূর্ণাঙ্ক (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, শূন্য *) = {
#ifdef XDOUBLE
  xtpmv_NUU, xtpmv_NUN, xtpmv_NLU, xtpmv_NLN,
  xtpmv_TUU, xtpmv_TUN, xtpmv_TLU, xtpmv_TLN,
  xtpmv_RUU, xtpmv_RUN, xtpmv_RLU, xtpmv_RLN,
  xtpmv_CUU, xtpmv_CUN, xtpmv_CLU, xtpmv_CLN,
#নোচেদ্যদি defined(DOUBLE)
  ztpmv_NUU, ztpmv_NUN, ztpmv_NLU, ztpmv_NLN,
  ztpmv_TUU, ztpmv_TUN, ztpmv_TLU, ztpmv_TLN,
  ztpmv_RUU, ztpmv_RUN, ztpmv_RLU, ztpmv_RLN,
  ztpmv_CUU, ztpmv_CUN, ztpmv_CLU, ztpmv_CLN,
#অন্যথা
  ctpmv_NUU, ctpmv_NUN, ctpmv_NLU, ctpmv_NLN,
  ctpmv_TUU, ctpmv_TUN, ctpmv_TLU, ctpmv_TLN,
  ctpmv_RUU, ctpmv_RUN, ctpmv_RLU, ctpmv_RLN,
  ctpmv_CUU, ctpmv_CUN, ctpmv_CLU, ctpmv_CLN,
#endif
};

#ifdef SMP
স্থির পূর্ণাঙ্ক (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, পূর্ণাঙ্ক) = {
#ifdef XDOUBLE
  xtpmv_thread_NUU, xtpmv_thread_NUN, xtpmv_thread_NLU, xtpmv_thread_NLN,
  xtpmv_thread_TUU, xtpmv_thread_TUN, xtpmv_thread_TLU, xtpmv_thread_TLN,
  xtpmv_thread_RUU, xtpmv_thread_RUN, xtpmv_thread_RLU, xtpmv_thread_RLN,
  xtpmv_thread_CUU, xtpmv_thread_CUN, xtpmv_thread_CLU, xtpmv_thread_CLN,
#নোচেদ্যদি defined(DOUBLE)
  ztpmv_thread_NUU, ztpmv_thread_NUN, ztpmv_thread_NLU, ztpmv_thread_NLN,
  ztpmv_thread_TUU, ztpmv_thread_TUN, ztpmv_thread_TLU, ztpmv_thread_TLN,
  ztpmv_thread_RUU, ztpmv_thread_RUN, ztpmv_thread_RLU, ztpmv_thread_RLN,
  ztpmv_thread_CUU, ztpmv_thread_CUN, ztpmv_thread_CLU, ztpmv_thread_CLN,
#অন্যথা
  ctpmv_thread_NUU, ctpmv_thread_NUN, ctpmv_thread_NLU, ctpmv_thread_NLN,
  ctpmv_thread_TUU, ctpmv_thread_TUN, ctpmv_thread_TLU, ctpmv_thread_TLN,
  ctpmv_thread_RUU, ctpmv_thread_RUN, ctpmv_thread_RLU, ctpmv_thread_RLN,
  ctpmv_thread_CUU, ctpmv_thread_CUN, ctpmv_thread_CLU, ctpmv_thread_CLN,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(বর্ণ *UPLO, বর্ণ *TRANS, বর্ণ *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  বর্ণ uplo_arg  = *UPLO;
  বর্ণ trans_arg = *TRANS;
  বর্ণ diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  পূর্ণাঙ্ক uplo;
  পূর্ণাঙ্ক unit;
  পূর্ণাঙ্ক trans;
  FLOAT *buffer;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'T') trans = 1;
  যদি (trans_arg == 'R') trans = 2;
  যদি (trans_arg == 'C') trans = 3;

  যদি (diag_arg  == 'U') unit  = 0;
  যদি (diag_arg  == 'N') unit  = 1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  যদি (incx == 0)          info =  7;
  যদি (n < 0)              info =  4;
  যদি (unit  < 0)          info =  3;
  যদি (trans < 0)          info =  2;
  যদি (uplo  < 0)          info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo,
	   গণনা CBLAS_TRANSPOSE TransA, গণনা CBLAS_DIAG Diag,
	   blasint n, শূন্য  *va, শূন্য  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  পূর্ণাঙ্ক trans, uplo, unit;
  blasint info;
  FLOAT *buffer;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;
#ifdef SMP
  পূর্ণাঙ্ক nthreads;
#endif

  যদি (order == CblasColMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 0;
    যদি (Uplo == CblasLower)         uplo  = 1;

    যদি (TransA == CblasNoTrans)     trans = 0;
    যদি (TransA == CblasTrans)       trans = 1;
    যদি (TransA == CblasConjNoTrans) trans = 2;
    যদি (TransA == CblasConjTrans)   trans = 3;

    যদি (Diag == CblasUnit)          unit  = 0;
    যদি (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    যদি (incx == 0)          info =  7;
    যদি (n < 0)              info =  4;
    যদি (unit  < 0)          info =  3;
    যদি (trans < 0)          info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (order == CblasRowMajor) {
    যদি (Uplo == CblasUpper)         uplo  = 1;
    যদি (Uplo == CblasLower)         uplo  = 0;

    যদি (TransA == CblasNoTrans)     trans = 1;
    যদি (TransA == CblasTrans)       trans = 0;
    যদি (TransA == CblasConjNoTrans) trans = 3;
    যদি (TransA == CblasConjTrans)   trans = 2;

    যদি (Diag == CblasUnit)          unit  = 0;
    যদি (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    যদি (incx == 0)          info =  7;
    যদি (n < 0)              info =  4;
    যদি (unit  < 0)          info =  3;
    যদি (trans < 0)          info =  2;
    যদি (uplo  < 0)          info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#endif

  যদি (n == 0) প্রত্যাগম;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#ifdef SMP
  nthreads = অঙ্ক_cpu_avail(2);

  যদি (n < 50) nthreads = 1;
  যদি (nthreads > 2 && n < 500) nthreads = 2;

  যদি (nthreads == 1) {
#endif

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer);

#ifdef SMP

  } অন্যথা {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  প্রত্যাগম;
}
