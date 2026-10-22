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
#परिभाषय ERROR_NAME "XSPR  "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZSPR  "
#अन्यथा
#परिभाषय ERROR_NAME "CSPR  "
#यदि_अन्तः

स्थिर पूर्णाङ्क (*spr[])(BLASLONG, FLOAT, FLOAT, FLOAT *, BLASLONG, FLOAT *, FLOAT *) = {
#यदि_परिभाषितम् XDOUBLE
  xspr_U, xspr_L,
#नोचेद्यदि defined(DOUBLE)
  zspr_U, zspr_L,
#अन्यथा
  cspr_U, cspr_L,
#यदि_अन्तः
};

#यदि_परिभाषितम् SMP
स्थिर पूर्णाङ्क (*spr_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, FLOAT *, पूर्णाङ्क) = {
#यदि_परिभाषितम् XDOUBLE
  xspr_thread_U, xspr_thread_L,
#नोचेद्यदि defined(DOUBLE)
  zspr_thread_U, zspr_thread_L,
#अन्यथा
  cspr_thread_U, cspr_thread_L,
#यदि_अन्तः
};
#यदि_अन्तः

OPENBLAS_EXPORT
शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT  *ALPHA,
	 FLOAT  *x, blasint *INCX, FLOAT *a){

  वर्ण uplo_arg = *UPLO;
  blasint n		= *N;
  FLOAT alpha_r	= ALPHA[0];
  FLOAT alpha_i	= ALPHA[1];
  blasint incx	= *INCX;

  blasint info;
  पूर्णाङ्क uplo;
  FLOAT *बफर;
#यदि_परिभाषितम् SMP
  पूर्णाङ्क nthreads;
#यदि_अन्तः

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  uplo  = -1;

  यदि (uplo_arg  == 'U') uplo  = 0;
  यदि (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  यदि (incx == 0)          info =  5;
  यदि (n < 0)              info =  2;
  यदि (uplo  < 0)          info =  1;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

  यदि (n == 0) प्रत्यागम;

  यदि ((alpha_r == ZERO) && (alpha_i == ZERO)) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  यदि (incx < 0 ) x -= (n - 1) * incx * COMPSIZE;

  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);
  यदि (!बफर) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
    प्रत्यागम;
  }

#यदि_परिभाषितम् SMP
  nthreads = अङ्कः_cpu_avail(2);

  यदि (n < 50) nthreads = 1;
  यदि (nthreads > 2 && n < 500) nthreads = 2;

  यदि (nthreads == 1) {
#यदि_अन्तः

    (spr[uplo])(n, alpha_r, alpha_i, x, incx, a, बफर);

#यदि_परिभाषितम् SMP
  } अन्यथा {

    (spr_thread[uplo])(n, ALPHA, x, incx, a, बफर, nthreads);

  }
#यदि_अन्तः

  blas_स्मृतिः_free(बफर);

  FUNCTION_PROFILE_अन्तः(4, n * n / 2 + n, n * n);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
