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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#ifdef XDOUBLE
#define ERROR_NAME "QGBMV "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DGBMV "
#अन्यथा
#define ERROR_NAME "SGBMV "
#endif

स्थिर शून्य (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#ifdef XDOUBLE
      qgbmv_n, qgbmv_t,
#किंवाजर defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#अन्यथा
      sgbmv_n, sgbmv_t,
#endif
};

#ifdef SMP
स्थिर पूर्णांक (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#ifdef XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#किंवाजर defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#अन्यथा
      sgbmv_thread_n, sgbmv_thread_t,
#endif
};
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

शून्य NAME(अक्षर *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  अक्षर trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *buffer;
#ifdef SMP
  पूर्णांक nthreads;
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

  जर (trans == 'N') i = 0;
  जर (trans == 'T') i = 1;
  जर (trans == 'R') i = 0;
  जर (trans == 'C') i = 1;

  जर (incy == 0)	 info = 13;
  जर (incx == 0)	 info = 10;
  जर (lda < kl + ku + 1) info = 8;
  जर (kl < 0)		 info = 5;
  जर (ku < 0)		 info = 4;
  जर (n < 0)		 info = 3;
  जर (m < 0)		 info = 2;
  जर (i < 0)		 info = 1;

  trans = i;

  जर (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order,
	   गणना CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *buffer;
  blasint lenx, leny, info, t;
  पूर्णांक trans;
#ifdef SMP
  पूर्णांक nthreads;
#endif

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  जर (order == CblasColMajor) {
    जर (TransA == CblasNoTrans)     trans = 0;
    जर (TransA == CblasTrans)       trans = 1;
    जर (TransA == CblasConjNoTrans) trans = 0;
    जर (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    जर (incy == 0)	 info = 13;
    जर (incx == 0)	 info = 10;
    जर (lda < kl + ku + 1) info = 8;
    जर (kl < 0)		 info = 5;
    जर (ku < 0)		 info = 4;
    जर (n < 0)		 info = 3;
    जर (m < 0)		 info = 2;
    जर (trans < 0)	 info = 1;
  }

  जर (order == CblasRowMajor) {
    जर (TransA == CblasNoTrans)     trans = 1;
    जर (TransA == CblasTrans)       trans = 0;
    जर (TransA == CblasConjNoTrans) trans = 1;
    जर (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    जर (incy == 0)	 info = 13;
    जर (incx == 0)	 info = 10;
    जर (lda < kl + ku + 1) info = 8;
    जर (kl < 0)		 info = 5;
    जर (ku < 0)		 info = 4;
    जर (n < 0)		 info = 3;
    जर (m < 0)		 info = 2;
    जर (trans < 0)	 info = 1;
  }

  जर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#endif

  जर ((m==0) || (n==0)) परत;

  lenx = n;
  leny = m;
  जर (trans) lenx = m;
  जर (trans) leny = n;

  जर (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  जर (alpha == ZERO) परत;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  जर (incx < 0) x -= (lenx-1)*incx;
  जर (incy < 0) y -= (leny-1)*incy;

  buffer = (FLOAT *)blas_memory_alloc(1);
  जर (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
    परत;
  }

#ifdef SMP
  जर (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  अन्यथा
    nthreads = अंक_cpu_avail(2);

  जर (nthreads == 1) {
#endif

  (gbmv[(पूर्णांक)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer);

#ifdef SMP
  } अन्यथा {

    (gbmv_thread[(पूर्णांक)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, buffer, nthreads);

  }
#endif

  blas_memory_free(buffer);

  FUNCTION_PROFILE_END(1, m * n / 2 + n, m * n);

  IDEBUG_END;

  परत;
}
