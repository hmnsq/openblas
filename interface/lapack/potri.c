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
#परिभाषय ERROR_NAME "QPOTRI"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DPOTRI"
#अन्यथा
#परिभाषय ERROR_NAME "SPOTRI"
#यदि_अन्तः

स्थिर blasint (*trtri_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UN_SINGLE, TRTRI_LN_SINGLE,
};

स्थिर blasint (*lauum_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  LAUUM_U_SINGLE, LAUUM_L_SINGLE,
};

#यदि_परिभाषितम् SMP
स्थिर blasint (*trtri_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  TRTRI_UN_PARALLEL, TRTRI_LN_PARALLEL,
};

स्थिर blasint (*lauum_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) ={
  LAUUM_U_PARALLEL, LAUUM_L_PARALLEL,
};
#यदि_अन्तः

OPENBLAS_EXPORT
शून्य NAME(वर्ण *UPLO, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t तर्काः;

  blasint uplo_arg = *UPLO;
  blasint uplo;
  blasint info;
  FLOAT *बफर;
#यदि_परिभाषितम् PPC440
  बाह्य
#यदि_अन्तः
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  तर्काः.n    = *N;
  तर्काः.a    = (शून्य *)a;
  तर्काः.lda  = *ldA;

  TOUPPER(uplo_arg);

  uplo = -1;
  यदि (uplo_arg == 'U') uplo = 0;
  यदि (uplo_arg == 'L') uplo = 1;

  info  = 0;
  यदि (तर्काः.lda  < MAX(1,तर्काः.n)) info = 4;
  यदि (तर्काः.n    < 0)             info = 2;
  यदि (uplo < 0)                  info = 1;

  यदि (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME) - 1);
    *Info = - info;
    प्रत्यागम;
  }

  *Info = 0;

  यदि (तर्काः.n == 0) प्रत्यागम;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#यदि_अपरिभाषितम् PPC440
  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);
#यदि_अन्तः

#यदि_परिभाषितम् SMP
  तर्काः.common = रिक्त;
  यदि (तर्काः.n < 180)
    तर्काः.nthreads = 1;
  अन्यथा
  तर्काः.nthreads = अङ्कः_cpu_avail(4);

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

  info = (trtri_single[uplo])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

  यदि (!info) {
    info = (lauum_single[uplo])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);
  }

  *Info = info;

#यदि_परिभाषितम् SMP
  } अन्यथा {
    info = (trtri_parallel[uplo])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

    यदि (!info) {
      info = (lauum_parallel[uplo])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);
    }

    *Info = info;
  }
#यदि_अन्तः

#यदि_अपरिभाषितम् PPC440
  blas_स्मृतिः_free(बफर);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, .5 * तर्काः.n * तर्काः.n,
		          तर्काः.n * (1./3. + तर्काः.n * ( 1./2. + तर्काः.n * 1./6.))
		       +  तर्काः.n * (1./3. + तर्काः.n * (-1./2. + तर्काः.n * 1./6.)));

  IDEBUG_अन्तः;

  प्रत्यागम;
}
