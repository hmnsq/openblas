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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGBMV "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGBMV "
#अन्यथा
#परिभाषय ERROR_NAME "SGBMV "
#यदि_अन्तः

स्थिर शून्य (*gbmv[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, शून्य *) = {
#यदि_परिभाषितम् XDOUBLE
      qgbmv_n, qgbmv_t,
#नोचेद्यदि defined(DOUBLE)
      dgbmv_n, dgbmv_t,
#अन्यथा
      sgbmv_n, sgbmv_t,
#यदि_अन्तः
};

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*gbmv_thread[])(BLASLONG, BLASLONG, BLASLONG, BLASLONG, FLOAT,
		      FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#यदि_परिभाषितम् XDOUBLE
      qgbmv_thread_n, qgbmv_thread_t,
#नोचेद्यदि defined(DOUBLE)
      dgbmv_thread_n, dgbmv_thread_t,
#अन्यथा
      sgbmv_thread_n, sgbmv_thread_t,
#यदि_अन्तः
};
#यदि_अन्तः

OPENBLAS_EXPORT
#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *TRANS, blasint *M, blasint *N,
	 blasint *KU, blasint *KL,
	 FLOAT *ALPHA, FLOAT *a, blasint *LDA,
	 FLOAT *x, blasint *INCX,
	 FLOAT *BETA, FLOAT *y, blasint *INCY){

  वर्ण trans = *TRANS;
  blasint m = *M;
  blasint n = *N;
  blasint ku = *KU;
  blasint kl = *KL;
  blasint lda = *LDA;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT *बफर;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  FLOAT alpha = *ALPHA;
  FLOAT beta  = *BETA;

  blasint info;
  blasint lenx, leny;
  blasint i;

  PRINT_DEBUG_NAME;

  TOUPPER(trans);

  info = 0;

  i = -1;

  यदि (trans == 'N') i = 0;
  यदि (trans == 'T') i = 1;
  यदि (trans == 'R') i = 0;
  यदि (trans == 'C') i = 1;

  यदि (incy == 0)	 info = 13;
  यदि (incx == 0)	 info = 10;
  यदि (lda < kl + ku + 1) info = 8;
  यदि (kl < 0)		 info = 5;
  यदि (ku < 0)		 info = 4;
  यदि (n < 0)		 info = 3;
  यदि (m < 0)		 info = 2;
  यदि (i < 0)		 info = 1;

  trans = i;

  यदि (info != 0){
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः,
	   गणना CBLAS_TRANSPOSE TransA,
	   blasint m, blasint n,
	   blasint ku, blasint kl,
	   FLOAT alpha,
	   FLOAT  *a, blasint lda,
	   FLOAT  *x, blasint incx,
	   FLOAT beta,
	   FLOAT  *y, blasint incy){

  FLOAT *बफर;
  blasint lenx, leny, info, t;
  पूर्णाङ्क trans;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  trans = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 0;
    यदि (TransA == CblasConjTrans)   trans = 1;

    info = -1;

    यदि (incy == 0)	 info = 13;
    यदि (incx == 0)	 info = 10;
    यदि (lda < kl + ku + 1) info = 8;
    यदि (kl < 0)		 info = 5;
    यदि (ku < 0)		 info = 4;
    यदि (n < 0)		 info = 3;
    यदि (m < 0)		 info = 2;
    यदि (trans < 0)	 info = 1;
  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 1;
    यदि (TransA == CblasConjTrans)   trans = 0;

    info = -1;

    t = n;
    n = m;
    m = t;

    t  = ku;
    ku = kl;
    kl = t;

    यदि (incy == 0)	 info = 13;
    यदि (incx == 0)	 info = 10;
    यदि (lda < kl + ku + 1) info = 8;
    यदि (kl < 0)		 info = 5;
    यदि (ku < 0)		 info = 4;
    यदि (n < 0)		 info = 3;
    यदि (m < 0)		 info = 2;
    यदि (trans < 0)	 info = 1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  यदि ((m==0) || (n==0)) प्रत्यागम;

  lenx = n;
  leny = m;
  यदि (trans) lenx = m;
  यदि (trans) leny = n;

  यदि (beta != ONE) SCAL_K(leny, 0, 0, beta, y, blasabs(incy), रिक्त, 0, रिक्त, 0);

  यदि (alpha == ZERO) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0) x -= (lenx-1)*incx;
  यदि (incy < 0) y -= (leny-1)*incy;

  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_परिभाषितम् SMP
  यदि (m * n < 250000 || kl+ku < 15 )
    nthreads = 1;
  अन्यथा
    nthreads = अङ्कः_cpu_avail(2);

  यदि (nthreads == 1) {
#यदि_अन्तः

  (gbmv[(पूर्णाङ्क)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, बफर);

#यदि_परिभाषितम् SMP
  } अन्यथा {

    (gbmv_thread[(पूर्णाङ्क)trans])(m, n, kl, ku, alpha, a, lda, x, incx, y, incy, बफर, nthreads);

  }
#यदि_अन्तः

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(1, m * n / 2 + n, m * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
