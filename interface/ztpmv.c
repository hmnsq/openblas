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

#समावेश <stdio.h>
#समावेश <ctype.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XTPMV "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZTPMV "
#अन्यथा
#घोषणा ERROR_NAME "CTPMV "
#पूर्ण_अगर

अटल पूर्णांक (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, व्योम *) = {
#अगर_घोषित XDOUBLE
  xtpmv_NUU, xtpmv_NUN, xtpmv_NLU, xtpmv_NLN,
  xtpmv_TUU, xtpmv_TUN, xtpmv_TLU, xtpmv_TLN,
  xtpmv_RUU, xtpmv_RUN, xtpmv_RLU, xtpmv_RLN,
  xtpmv_CUU, xtpmv_CUN, xtpmv_CLU, xtpmv_CLN,
#या_अगर defined(DOUBLE)
  ztpmv_NUU, ztpmv_NUN, ztpmv_NLU, ztpmv_NLN,
  ztpmv_TUU, ztpmv_TUN, ztpmv_TLU, ztpmv_TLN,
  ztpmv_RUU, ztpmv_RUN, ztpmv_RLU, ztpmv_RLN,
  ztpmv_CUU, ztpmv_CUN, ztpmv_CLU, ztpmv_CLN,
#अन्यथा
  ctpmv_NUU, ctpmv_NUN, ctpmv_NLU, ctpmv_NLN,
  ctpmv_TUU, ctpmv_TUN, ctpmv_TLU, ctpmv_TLN,
  ctpmv_RUU, ctpmv_RUN, ctpmv_RLU, ctpmv_RLN,
  ctpmv_CUU, ctpmv_CUN, ctpmv_CLU, ctpmv_CLN,
#पूर्ण_अगर
};

#अगर_घोषित SMP
अटल पूर्णांक (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, पूर्णांक) = {
#अगर_घोषित XDOUBLE
  xtpmv_thread_NUU, xtpmv_thread_NUN, xtpmv_thread_NLU, xtpmv_thread_NLN,
  xtpmv_thread_TUU, xtpmv_thread_TUN, xtpmv_thread_TLU, xtpmv_thread_TLN,
  xtpmv_thread_RUU, xtpmv_thread_RUN, xtpmv_thread_RLU, xtpmv_thread_RLN,
  xtpmv_thread_CUU, xtpmv_thread_CUN, xtpmv_thread_CLU, xtpmv_thread_CLN,
#या_अगर defined(DOUBLE)
  ztpmv_thread_NUU, ztpmv_thread_NUN, ztpmv_thread_NLU, ztpmv_thread_NLN,
  ztpmv_thread_TUU, ztpmv_thread_TUN, ztpmv_thread_TLU, ztpmv_thread_TLN,
  ztpmv_thread_RUU, ztpmv_thread_RUN, ztpmv_thread_RLU, ztpmv_thread_RLN,
  ztpmv_thread_CUU, ztpmv_thread_CUN, ztpmv_thread_CLU, ztpmv_thread_CLN,
#अन्यथा
  ctpmv_thread_NUU, ctpmv_thread_NUN, ctpmv_thread_NLU, ctpmv_thread_NLN,
  ctpmv_thread_TUU, ctpmv_thread_TUN, ctpmv_thread_TLU, ctpmv_thread_TLN,
  ctpmv_thread_RUU, ctpmv_thread_RUN, ctpmv_thread_RLU, ctpmv_thread_RLN,
  ctpmv_thread_CUU, ctpmv_thread_CUN, ctpmv_thread_CLU, ctpmv_thread_CLN,
#पूर्ण_अगर
};
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(अक्षर *UPLO, अक्षर *TRANS, अक्षर *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  अक्षर uplo_arg  = *UPLO;
  अक्षर trans_arg = *TRANS;
  अक्षर diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  पूर्णांक uplo;
  पूर्णांक unit;
  पूर्णांक trans;
  FLOAT *बफर;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  अगर (trans_arg == 'N') trans = 0;
  अगर (trans_arg == 'T') trans = 1;
  अगर (trans_arg == 'R') trans = 2;
  अगर (trans_arg == 'C') trans = 3;

  अगर (diag_arg  == 'U') unit  = 0;
  अगर (diag_arg  == 'N') unit  = 1;

  अगर (uplo_arg  == 'U') uplo  = 0;
  अगर (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  अगर (incx == 0)          info =  7;
  अगर (n < 0)              info =  4;
  अगर (unit  < 0)          info =  3;
  अगर (trans < 0)          info =  2;
  अगर (uplo  < 0)          info =  1;

  अगर (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_UPLO Uplo,
	   क्रमागत CBLAS_TRANSPOSE TransA, क्रमागत CBLAS_DIAG Diag,
	   blasint n, व्योम  *va, व्योम  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
  
  पूर्णांक trans, uplo, unit;
  blasint info;
  FLOAT *बफर;

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;
#अगर_घोषित SMP
  पूर्णांक nthreads;
#पूर्ण_अगर

  अगर (order == CblasColMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 0;
    अगर (Uplo == CblasLower)         uplo  = 1;

    अगर (TransA == CblasNoTrans)     trans = 0;
    अगर (TransA == CblasTrans)       trans = 1;
    अगर (TransA == CblasConjNoTrans) trans = 2;
    अगर (TransA == CblasConjTrans)   trans = 3;

    अगर (Diag == CblasUnit)          unit  = 0;
    अगर (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    अगर (incx == 0)          info =  7;
    अगर (n < 0)              info =  4;
    अगर (unit  < 0)          info =  3;
    अगर (trans < 0)          info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (order == CblasRowMajor) {
    अगर (Uplo == CblasUpper)         uplo  = 1;
    अगर (Uplo == CblasLower)         uplo  = 0;

    अगर (TransA == CblasNoTrans)     trans = 1;
    अगर (TransA == CblasTrans)       trans = 0;
    अगर (TransA == CblasConjNoTrans) trans = 3;
    अगर (TransA == CblasConjTrans)   trans = 2;

    अगर (Diag == CblasUnit)          unit  = 0;
    अगर (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    अगर (incx == 0)          info =  7;
    अगर (n < 0)              info =  4;
    अगर (unit  < 0)          info =  3;
    अगर (trans < 0)          info =  2;
    अगर (uplo  < 0)          info =  1;
  }

  अगर (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#पूर्ण_अगर

  अगर (n == 0) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  अगर (incx < 0 ) x -= (n - 1) * incx * 2;

  बफर = (FLOAT *)blas_memory_alloc(1);
  अगर (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
    वापस;
  }

#अगर_घोषित SMP
  nthreads = अंक_cpu_avail(2);

  अगर (n < 50) nthreads = 1;
  अगर (nthreads > 2 && n < 500) nthreads = 2;

  अगर (nthreads == 1) {
#पूर्ण_अगर

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, बफर);

#अगर_घोषित SMP

  } अन्यथा {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, बफर, nthreads);

  }
#पूर्ण_अगर

  blas_memory_free(बफर);

  FUNCTION_PROFILE_END(4, n * n / 2 + n, n * n);

  IDEBUG_END;

  वापस;
}
