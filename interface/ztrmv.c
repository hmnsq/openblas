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
#अन्तर्भाव <ctype.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XTRMV "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZTRMV "
#अन्यथा
#परिभाषय ERROR_NAME "CTRMV "
#यदि_अन्तः

स्थिर पूर्णाङ्क (*trmv[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *) = {
#यदि_परिभाषितम् XDOUBLE
  xtrmv_NUU, xtrmv_NUN, xtrmv_NLU, xtrmv_NLN,
  xtrmv_TUU, xtrmv_TUN, xtrmv_TLU, xtrmv_TLN,
  xtrmv_RUU, xtrmv_सततम्, xtrmv_RLU, xtrmv_RLN,
  xtrmv_CUU, xtrmv_CUN, xtrmv_CLU, xtrmv_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztrmv_NUU, ztrmv_NUN, ztrmv_NLU, ztrmv_NLN,
  ztrmv_TUU, ztrmv_TUN, ztrmv_TLU, ztrmv_TLN,
  ztrmv_RUU, ztrmv_सततम्, ztrmv_RLU, ztrmv_RLN,
  ztrmv_CUU, ztrmv_CUN, ztrmv_CLU, ztrmv_CLN,
#अन्यथा
  ctrmv_NUU, ctrmv_NUN, ctrmv_NLU, ctrmv_NLN,
  ctrmv_TUU, ctrmv_TUN, ctrmv_TLU, ctrmv_TLN,
  ctrmv_RUU, ctrmv_सततम्, ctrmv_RLU, ctrmv_RLN,
  ctrmv_CUU, ctrmv_CUN, ctrmv_CLU, ctrmv_CLN,
#यदि_अन्तः
};

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*trmv_thread[])(BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, पूर्णाङ्क) = {
#यदि_परिभाषितम् XDOUBLE
  xtrmv_thread_NUU, xtrmv_thread_NUN, xtrmv_thread_NLU, xtrmv_thread_NLN,
  xtrmv_thread_TUU, xtrmv_thread_TUN, xtrmv_thread_TLU, xtrmv_thread_TLN,
  xtrmv_thread_RUU, xtrmv_thread_सततम्, xtrmv_thread_RLU, xtrmv_thread_RLN,
  xtrmv_thread_CUU, xtrmv_thread_CUN, xtrmv_thread_CLU, xtrmv_thread_CLN,
#नोचेद्यदि defined(DOUBLE)
  ztrmv_thread_NUU, ztrmv_thread_NUN, ztrmv_thread_NLU, ztrmv_thread_NLN,
  ztrmv_thread_TUU, ztrmv_thread_TUN, ztrmv_thread_TLU, ztrmv_thread_TLN,
  ztrmv_thread_RUU, ztrmv_thread_सततम्, ztrmv_thread_RLU, ztrmv_thread_RLN,
  ztrmv_thread_CUU, ztrmv_thread_CUN, ztrmv_thread_CLU, ztrmv_thread_CLN,
#अन्यथा
  ctrmv_thread_NUU, ctrmv_thread_NUN, ctrmv_thread_NLU, ctrmv_thread_NLN,
  ctrmv_thread_TUU, ctrmv_thread_TUN, ctrmv_thread_TLU, ctrmv_thread_TLN,
  ctrmv_thread_RUU, ctrmv_thread_सततम्, ctrmv_thread_RLU, ctrmv_thread_RLN,
  ctrmv_thread_CUU, ctrmv_thread_CUN, ctrmv_thread_CLU, ctrmv_thread_CLN,
#यदि_अन्तः
};
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(वर्ण *UPLO, वर्ण *TRANS, वर्ण *DIAG,
	   blasint *N, FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  वर्ण uplo_arg  = *UPLO;
  वर्ण trans_arg = *TRANS;
  वर्ण diag_arg  = *DIAG;

  blasint n    = *N;
  blasint lda  = *LDA;
  blasint incx = *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  पूर्णाङ्क unit;
  पूर्णाङ्क trans, बफर_आकार;
  FLOAT *बफर;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'R') trans = 2;
  यदि (trans_arg == 'C') trans = 3;

  यदि (diag_arg  == 'U') unit  = 0;
  यदि (diag_arg  == 'N') unit  = 1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  यदि (incx == 0)          info =  8;
  यदि (lda  < MAX(1, n))   info =  6;
  यदि (n < 0)              info =  4;
  यदि (unit  < 0)          info =  3;
  यदि (trans < 0)          info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_UPLO Uplo,
	   गणना CBLAS_TRANSPOSE TransA, गणना CBLAS_DIAG Diag,
	   blasint n, शून्य  *va, blasint lda, शून्य  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  पूर्णाङ्क trans, uplo, unit, बफर_आकार;
  blasint info;
  FLOAT *बफर;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  यदि (क्रमः == CblasColMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)         uplo  = 0;
    यदि (Uplo == Cblas_अधस्तनम्)         uplo  = 1;

    यदि (TransA == CblasNoTrans)     trans = 0;
    यदि (TransA == CblasTrans)       trans = 1;
    यदि (TransA == CblasConjNoTrans) trans = 2;
    यदि (TransA == CblasConjTrans)   trans = 3;

    यदि (Diag == CblasUnit)          unit  = 0;
    यदि (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    यदि (incx == 0)          info =  8;
    यदि (lda  < MAX(1, n))   info =  6;
    यदि (n < 0)              info =  4;
    यदि (unit  < 0)          info =  3;
    यदि (trans < 0)          info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (क्रमः == CblasRowMajor) {
    यदि (Uplo == Cblas_ऊर्ध्वस्थम्)         uplo  = 1;
    यदि (Uplo == Cblas_अधस्तनम्)         uplo  = 0;

    यदि (TransA == CblasNoTrans)     trans = 1;
    यदि (TransA == CblasTrans)       trans = 0;
    यदि (TransA == CblasConjNoTrans) trans = 3;
    यदि (TransA == CblasConjTrans)   trans = 2;

    यदि (Diag == CblasUnit)          unit  = 0;
    यदि (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    यदि (incx == 0)          info =  8;
    यदि (lda  < MAX(1, n))   info =  6;
    यदि (n < 0)              info =  4;
    यदि (unit  < 0)          info =  3;
    यदि (trans < 0)          info =  2;
    यदि (uplo  < 0)          info =  1;
  }

  यदि (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_अन्तः

  यदि (n == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx * 2;

#यदि_परिभाषितम् SMP
  // Calibrated on a Xeon E5-2630
  यदि(1L * n * n > 36L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_देहली) {
    nthreads = अङ्कः_cpu_avail(2);
    यदि(nthreads > 2 && 1L * n * n < 64L * परिमाण(FLOAT) * परिमाण(FLOAT) * GEMM_MULTITHREAD_देहली)
      nthreads = 2;
  } अन्यथा
      nthreads = 1;

  यदि(nthreads > 1) {
    बफर_आकार = n > 16 ? 0 : n * 4 + 40;
  }
  अन्यथा
#यदि_अन्तः
  {
    बफर_आकार = ((n - 1) / DTB_ENTRIES) * 2 * DTB_ENTRIES + 32 / परिमाण(FLOAT);
    // It seems to be required for some K8 or Barcelona CPU
    बफर_आकार += 8;
    यदि(incx != 1)
      बफर_आकार += n * 2;
  }
  राशिः_ALLOC(बफर_आकार, FLOAT, बफर);

#यदि_परिभाषितम् SMP
  यदि (nthreads == 1) {
#यदि_अन्तः

  (trmv[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, बफर);

#यदि_परिभाषितम् SMP
  } अन्यथा {

    (trmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, lda, x, incx, बफर, nthreads);

  }
#यदि_अन्तः

  राशिः_FREE(बफर);

  FUNCTION_PROFILE_अन्तः(4, n * n / 2 + n,  n * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
