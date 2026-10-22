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
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#elif defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#မဟုတ်ပါက
#define ERROR_NAME "SGBMV "
#endif

static void (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, void *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#elif defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#မဟုတ်ပါက
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
static ကိန်းပြည့် (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ကိန်းပြည့်) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#elif defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#မဟုတ်ပါက
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

void NAME(var_na *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  var_na trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  ကိန်းပြည့် nthreads;
#endif

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  အကယ်၍ (trans == 'N') i = 0;
  အကယ်၍ (trans == 'T') i = 1;
  အကယ်၍ (trans == 'R') i = 0;
  အကယ်၍ (trans == 'C') i = 1;

  အကယ်၍ (incy == 0)	 info = 13;
  အကယ်၍ (incx == 0)	 info = 10;
  အကယ်၍ (lda < kl + ku + 1) info = 8;
  အကယ်၍ (kl < 0)		 info = 5;
  အကယ်၍ (ku < 0)		 info = 4;
  အကယ်၍ (n < 0)		 info = 3;
  အကယ်၍ (m < 0)		 info = 2;
  အကယ်၍ (i < 0)		 info = 1;

  trans = i;

  အကယ်၍ (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#မဟုတ်ပါက

void CNAME(ga_nanaa CBLAS_ORDER order,
	   ga_nanaa CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  ကိန်းပြည့် trans;
#ifdef SMP
  ကိန်းပြည့် nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  အကယ်၍ (order == CblasColMajor) {
    အကယ်၍ (TransA == CblasNoTrans)     trans = 0;
    အကယ်၍ (TransA == CblasTrans)       trans = 1;
    အကယ်၍ (TransA == CblasConjNoTrans) trans = 0;
    အကယ်၍ (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    အကယ်၍ (incy == 0)	 info = 13;
    အကယ်၍ (incx == 0)	 info = 10;
    အကယ်၍ (lda < kl + ku + 1) info = 8;
    အကယ်၍ (kl < 0)		 info = 5;
    အကယ်၍ (ku < 0)		 info = 4;
    အကယ်၍ (n < 0)		 info = 3;
    အကယ်၍ (m < 0)		 info = 2;
    အကယ်၍ (trans < 0)	 info = 1;
  }

  အကယ်၍ (order == CblasRowMajor) {
    အကယ်၍ (TransA == CblasNoTrans)     trans = 1;
    အကယ်၍ (TransA == CblasTrans)       trans = 0;
    အကယ်၍ (TransA == CblasConjNoTrans) trans = 1;
    အကယ်၍ (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    အကယ်၍ (incy == 0)	 info = 13;
    အကယ်၍ (incx == 0)	 info = 10;
    အကယ်၍ (lda < kl + ku + 1) info = 8;
    အကယ်၍ (kl < 0)		 info = 5;
    အကယ်၍ (ku < 0)		 info = 4;
    အကယ်၍ (n < 0)		 info = 3;
    အကယ်၍ (m < 0)		 info = 2;
    အကယ်၍ (trans < 0)	 info = 1;
  }

  အကယ်၍ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#endif

  အကယ်၍ ((m==0) || (n==0)) ပြန်ပေး;

  lenx = n;
  leny = m;
  အကယ်၍ (trans) lenx = m;
  အကယ်၍ (trans) leny = n;

  အကယ်၍ (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), NULL, 0, NULL, 0);

  အကယ်၍ (alpha == ZERO) ပြန်ပေး;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  အကယ်၍ (incx < 0) x -= (lenx-1)*incx;
  အကယ်၍ (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  အကယ်၍ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#ifdef SMP
  အကယ်၍ (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  မဟုတ်ပါက
    nthreads = num_cpu_avail(2);

  အကယ်၍ (nthreads == 1) {
#endif

  (gbmv[(ကိန်းပြည့်)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } မဟုတ်ပါက {

    (gbmv_thread[(ကိန်းပြည့်)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  ပြန်ပေး;
}
