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
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYR  "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYR  "
#ورنہ
#تعریف_کرو ERROR_NAME "CSYR  "
#اگر_ختم

ساکن صحیح (*syr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#اگر_تعریف_شدہ XDOUBLE
  xsyr_U, xsyr_L,
#وگرنہ_اگر defined(DOUBLE)
  zsyr_U, zsyr_L,
#ورنہ
  csyr_U, csyr_L,
#اگر_ختم
};

#اگر_تعریف_شدہ SMP
ساکن صحیح (*syr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#اگر_تعریف_شدہ XDOUBLE
  xsyr_thread_U, xsyr_thread_L,
#وگرنہ_اگر defined(DOUBLE)
  zsyr_thread_U, zsyr_thread_L,
#ورنہ
  csyr_thread_U, csyr_thread_L,
#اگر_ختم
};
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a, blasint *LDA){

  حرف uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint lda	= *LDA;
  blasint incx	= *INCX;

  blasint info;
  صحیح uplo;
  FLOAT *بفر;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  اگر (lda  < MAX(1, n))   info =  7;
  اگر (incx == 0)          info =  5;
  اگر (n < 0)              info =  2;
  اگر (uplo  < 0)          info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }


#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo, صحیح n, خالی* valpha, FLOAT *x, صحیح incx, FLOAT *a, صحیح lda) {

  FLOAT *بفر;
  صحیح uplo;
  blasint info;
  FLOAT * ALPHA = (FLOAT*)valpha;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_CNAME;

  uplo  = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {

    اگر (Uplo == Cblas_بالائی) uplo  = 0;
    اگر (Uplo == Cblas_زیریں) uplo  = 1;

    info = -1;

    اگر (lda  < MAX(1, n))   info =  7;
    اگر (incx == 0)          info =  5;
    اگر (n < 0)              info =  2;
    اگر (uplo  < 0)          info =  1;

  }

  اگر (ترتیب == Cblas_صف_Major) {

    اگر (Uplo == Cblas_بالائی) uplo  = 1;
    اگر (Uplo == Cblas_زیریں) uplo  = 0;

    info = -1;

    اگر (lda  < MAX(1, n))   info =  7;
    اگر (incx == 0)          info =  5;
    اگر (n < 0)              info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_ختم

  اگر (n == 0) واپس;

  اگر ((alpha_r == ZERO) && (alpha_i == ZERO)) واپس;

  اگر (incx == 1 && n < 50) {
    blasint i;
    اگر (!uplo) {
      برائے (i = 0; i < n; i++){
        اگر ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(i + 1, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x,         1, a, 1, تہی, 0);
        }
        a += lda * 2;
      }
    } ورنہ {
      برائے (i = 0; i < n; i++){
        اگر ((x[i * 2 + 0] != ZERO) || (x[i * 2 + 1] != ZERO)) {
          AXPYU_K(n - i, 0, 0,
              alpha_r * x[i * 2 + 0] - alpha_i * x[i * 2 + 1],
              alpha_i * x[i * 2 + 0] + alpha_r * x[i * 2 + 1],
              x + i * 2, 1, a, 1, تہی, 0);
        }
        a += 2 + lda * 2;
      }
    }
    واپس;
  }

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0 ) x -= (n - 1) * incx * 2;

  بفر = (FLOAT *)blas_حافظہ_alloc(1);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_تعریف_شدہ SMP
  nthreads = عدد_cpu_avail(2);

  اگر (n < 50) nthreads = 1;
  اگر (nthreads > 2 && n < 500) nthreads = 2;

  اگر (nthreads == 1) {
#اگر_ختم

  (syr[uplo])(n, alpha_r, alpha_i, x, incx, a, lda, بفر);

#اگر_تعریف_شدہ SMP
  } ورنہ {

    (syr_thread[uplo])(n, ALPHA, x, incx, a, lda, بفر, nthreads);

  }
#اگر_ختم

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(4, n * n / 2 + n, n * n);

  IDEBUG_انتہا;

  واپس;
}
