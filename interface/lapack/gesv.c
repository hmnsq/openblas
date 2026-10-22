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

#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGESV"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGESV"
#अन्यथा
#परिभाषय ERROR_NAME "SGESV"
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGESV"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGESV"
#अन्यथा
#परिभाषय ERROR_NAME "CGESV"
#यदि_अन्तः
#यदि_अन्तः

OPENBLAS_EXPORT
शून्य NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t तर्काः;

  blasint info;
  FLOAT *बफर;
#यदि_परिभाषितम् PPC440
  बाह्य
#यदि_अन्तः
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  तर्काः.m    = *N;
  तर्काः.n    = *NRHS;
  तर्काः.a    = (शून्य *)a;
  तर्काः.lda  = *ldA;
  तर्काः.b    = (शून्य *)b;
  तर्काः.ldb  = *ldB;
  तर्काः.c    = (शून्य *)ipiv;

  info  = 0;
  यदि (तर्काः.ldb < MAX(1,तर्काः.m)) info = 7;
  यदि (तर्काः.lda < MAX(1,तर्काः.m)) info = 4;
  यदि (तर्काः.n   < 0)             info = 2;
  यदि (तर्काः.m   < 0)             info = 1;

  यदि (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME) - 1);
    *Info = - info;
    प्रत्यागम;
  }

  तर्काः.alpha = रिक्त;
  तर्काः.beta  = रिक्त;

  *Info = 0;

  यदि (तर्काः.m == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#यदि_अपरिभाषितम् PPC440
  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);
#यदि_अन्तः

#यदि_परिभाषितम् SMP
  तर्काः.common = रिक्त;

#यदि defined(_WIN64) && defined(_M_ARM64)
  #यदि_परिभाषितम् COMPLEX
    यदि (तर्काः.m * तर्काः.m <= 300) 
  #अन्यथा
    यदि (तर्काः.m * तर्काः.m <= 500) 
  #यदि_अन्तः
      तर्काः.nthreads = 1;
  अन्यथा यदि (तर्काः.m * तर्काः.m <= 1000)
    तर्काः.nthreads = 4;
  अन्यथा
    तर्काः.nthreads = अङ्कः_cpu_avail(4);
#अन्यथा
  #यदि_अपरिभाषितम् DOUBLE
    यदि (तर्काः.m * तर्काः.m < 40000)
  #अन्यथा
    यदि (तर्काः.m * तर्काः.m < 10000)
  #यदि_अन्तः
      तर्काः.nthreads = 1;
    अन्यथा
      तर्काः.nthreads = अङ्कः_cpu_avail(4);
#यदि_अन्तः

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

    तर्काः.n    = *N;
    info = GETRF_SINGLE(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

    यदि (info == 0 && *NRHS >0){
      तर्काः.n    = *NRHS;
      GETRS_N_SINGLE(&तर्काः, रिक्त, रिक्त, sa, sb, 0);
    }

#यदि_परिभाषितम् SMP
  } अन्यथा {

    तर्काः.n    = *N;
    info = GETRF_PARALLEL(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

    यदि (info == 0 && *NRHS > 0){
      तर्काः.n    = *NRHS;
      GETRS_N_PARALLEL(&तर्काः, रिक्त, रिक्त, sa, sb, 0);
    }
  }
#यदि_अन्तः

#यदि_अपरिभाषितम् PPC440
  blas_स्मृतिः_free(बफर);
#यदि_अन्तः

  *Info = info;

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_अन्तः;

  प्रत्यागम;
}
