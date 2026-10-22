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
#define ERROR_NAME "QSPR  "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DSPR  "
#નહીંતર
#define ERROR_NAME "SSPR  "
#endif

સ્થિર પૂર્ણાંક (*spr[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#ifdef XDOUBLE
  qspr_U, qspr_L,
#અથવાજો defined(DOUBLE)
  dspr_U, dspr_L,
#નહીંતર
  sspr_U, sspr_L,
#endif
};

#ifdef SMP
સ્થિર પૂર્ણાંક (*spr_thread[])(BLASLONG, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *, પૂર્ણાંક) = {
#ifdef XDOUBLE
  qspr_thread_U, qspr_thread_L,
#અથવાજો defined(DOUBLE)
  dspr_thread_U, dspr_thread_L,
#નહીંતર
  sspr_thread_U, sspr_thread_L,
#endif
};
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  અક્ષર uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha	= *ALPHA;
  blasint incx	= *INCX;

  blasint info;
  પૂર્ણાંક uplo;
  FLOAT *buffer;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  જો (incx == 0)          info =  5;
  જો (n < 0)              info =  2;
  જો (uplo  < 0)          info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order,
	   ગણના CBLAS_UPLO Uplo,
	   blasint n,
	   FLOAT alpha,
	   FLOAT  *x, blasint incx,
	   FLOAT  *a) {

  FLOAT *buffer;
  પૂર્ણાંક uplo;
  blasint info;
#ifdef SMP
  પૂર્ણાંક nthreads;
#endif

  PRINT_DEBUG_CNAME;

  uplo = -1;
  info =  0;

  જો (order == CblasColMajor) {
    જો (Uplo == CblasUpper)         uplo  = 0;
    જો (Uplo == CblasLower)         uplo  = 1;

    info = -1;

    જો (incx == 0)          info =  5;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Uplo == CblasUpper)         uplo  = 1;
    જો (Uplo == CblasLower)         uplo  = 0;

    info = -1;

    જો (incx == 0)          info =  5;
    જો (n < 0)              info =  2;
    જો (uplo  < 0)          info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#endif

  જો (n == 0) પાછા;

  જો (alpha == ZERO) પાછા;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  જો (incx == 1 && n <100) {
    blasint i;
    જો (uplo==0) {
      માટે (i = 0; i < n; i++){
        જો (x[i] != ZERO) {
          AXPYU_K(i + 1, 0, 0, alpha * x[i], x,     1, a, 1, ખાલી, 0);
        }
        a += i + 1;
      }
    } નહીંતર { 
      માટે (i = 0; i < n; i++){
        જો (x[i] != ZERO) {
          AXPYU_K(n - i, 0, 0, alpha * x[i], x + i, 1, a, 1, ખાલી, 0);
        }
        a += n - i;
      }
    }
    પાછા;
  }

  જો (incx < 0 ) x -= (n - 1) * incx;

  buffer = (FLOAT *)blas_memory_alloc(1);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#ifdef SMP
  nthreads = અંક_cpu_avail(2);

  જો (n < 50) nthreads = 1;
  જો (nthreads > 2 && n < 500) nthreads = 2;

  જો (nthreads == 1) {
#endif

  (spr[uplo])(n, alpha, x, incx, a, buffer);

#ifdef SMP
  } નહીંતર {

    (spr_thread[uplo])(n, alpha, x, incx, a, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, n * n / 2 + n, n * n);

  IDEBUG_END;

  પાછા;
}
