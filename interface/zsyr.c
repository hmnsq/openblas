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
#define ERROR_NAME "XSYR  "
#elif defined(DOUBLE)
#define ERROR_NAME "ZSYR  "
#မဟုတ်ပါက
#define ERROR_NAME "CSYR  "
#endif

static ကိန်းပြည့် (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#ifdef XDOUBLE
  xsyr_U, xsyr_L,
#elif defined(DOUBLE)
  zsyr_U, zsyr_L,
#မဟုတ်ပါက
  csyr_U, csyr_L,
#endif
};

#ifdef SMP
static ကိန်းပြည့် (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, ကိန်းပြည့်) = {
#ifdef XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#elif defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#မဟုတ်ပါက
  csyr_thread_U, csyr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(var_na *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  var_na uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  ကိန်းပြည့် uplo;
  FLOAT *buffer;
#ifdef SMP
  ကိန်းပြည့် nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  အကယ်၍ (uplo_arg  == 'U') uplo  = 0;
  အကယ်၍ (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  အကယ်၍ (lda  < MAX(1, n))   info =  7;
  အကယ်၍ (incx == 0)          info =  5;
  အကယ်၍ (n < 0)              info =  2;
  အကယ်၍ (uplo  < 0)          info =  1;

  အကယ်၍ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }


#မဟုတ်ပါက

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_UPLO Uplo, ကိန်းပြည့် n, void* valpha, FLOAT *x, ကိန်းပြည့် incx, FLOAT *a, ကိန်းပြည့် lda) {

  FLOAT *buffer;
  ကိန်းပြည့် uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
#ifdef SMP
  ကိန်းပြည့် nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  အကယ်၍ (order == CblasColMajor) {

    အကယ်၍ (Uplo == CblasUpper) uplo  = 0;
    အကယ်၍ (Uplo == CblasLower) uplo  = 1;

    info = -1;

    အကယ်၍ (lda  < MAX(1, n))   info =  7;
    အကယ်၍ (incx == 0)          info =  5;
    အကယ်၍ (n < 0)              info =  2;
    အကယ်၍ (uplo  < 0)          info =  1;

  }

  အကယ်၍ (order == CblasRowMajor) {

    အကယ်၍ (Uplo == CblasUpper) uplo  = 1;
    အကယ်၍ (Uplo == CblasLower) uplo  = 0;

    info = -1;

    အကယ်၍ (lda  < MAX(1, n))   info =  7;
    အကယ်၍ (incx == 0)          info =  5;
    အကယ်၍ (n < 0)              info =  2;
    အကယ်၍ (uplo  < 0)          info =  1;
  }

  အကယ်၍ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#endif

  အကယ်၍ (n == 0) ပြန်ပေး;

  အကယ်၍ ((alpha_r == ZERO) && (alpha_i == ZERO)) ပြန်ပေး;

  အကယ်၍ (incx == 1 && n < 50) {
    blasint i;
    အကယ်၍ (!uplo) {
      အတွက် (i = 0; i < n; i++){
        အကယ်၍ ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, NULL, 0);
        }
        a += lda * 2;
      }
    } မဟုတ်ပါက {
      အတွက် (i = 0; i < n; i++){
        အကယ်၍ ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, NULL, 0);
        }
        a += 2 + lda * 2;
      }
    }
    ပြန်ပေး;
  }

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  အကယ်၍ (incx < 0 ) x -= (n - 1) * incx * 2;

  buffer = (FLOAT *)blas_memory_alloc(1);
  အကယ်၍ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
    ပြန်ပေး;
  }

#ifdef SMP
  nthreads = num_cpu_avail(2);

  အကယ်၍ (n < 50) nthreads = 1;
  အကယ်၍ (nthreads > 2 && n < 500) nthreads = 2;

  အကယ်၍ (nthreads == 1) {
#endif

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, buffer);

#ifdef SMP
  } မဟုတ်ပါက {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  ပြန်ပေး;
}
