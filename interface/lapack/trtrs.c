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
#परिभाषय ERROR_NAME "QTRTRS"
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DTRTRS"
#अन्यथा
#परिभाषय ERROR_NAME "STRTRS"
#यदि_अन्तः

स्थिर blasint (*trtrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
    TRTRS_UNU_SINGLE, TRTRS_UNN_SINGLE, TRTRS_UTU_SINGLE, TRTRS_UTN_SINGLE, TRTRS_LNU_SINGLE, TRTRS_LNN_SINGLE, TRTRS_LTU_SINGLE, TRTRS_LTN_SINGLE,
};

#यदि_परिभाषितम् SMP
स्थिर blasint (*trtrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
    TRTRS_UNU_PARALLEL, TRTRS_UNN_PARALLEL, TRTRS_UTU_PARALLEL, TRTRS_UTN_PARALLEL, TRTRS_LNU_PARALLEL, TRTRS_LNN_PARALLEL, TRTRS_LTU_PARALLEL, TRTRS_LTN_PARALLEL,
};
#यदि_अन्तः

OPENBLAS_EXPORT
शून्य NAME(वर्ण *UPLO, वर्ण* TRANS, वर्ण* DIAG, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  FLOAT *b, blasint *ldB, blasint *Info){

    वर्ण uplo_arg = *UPLO;
    वर्ण trans_arg = *TRANS;
    वर्ण diag_arg = *DIAG;

  blas_arg_t तर्काः;

  blasint info;
  पूर्णाङ्क uplo, trans, diag;
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

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;
  यदि (trans_arg == 'N') trans = 0;
  यदि (trans_arg == 'T') trans = 1;
  यदि (trans_arg == 'R') trans = 0;
  यदि (trans_arg == 'C') trans = 1;

  TOUPPER(uplo_arg);
  uplo = -1;
  यदि (uplo_arg == 'U') uplo = 0;
  यदि (uplo_arg == 'L') uplo = 1;

  TOUPPER(diag_arg);
  diag = -1;
  यदि (diag_arg == 'U') diag = 0;
  यदि (diag_arg == 'N') diag = 1;

  यदि (तर्काः.ldb  < MAX(1, तर्काः.m)) info = 9;
  यदि (तर्काः.lda  < MAX(1, तर्काः.m)) info = 7;
  यदि (तर्काः.n    < 0) info = 5;
  यदि (तर्काः.m    < 0) info = 4;
  यदि (trans     < 0) info = 2;
  यदि (uplo      < 0) info = 1;
  यदि (diag      < 0) info = 3;

  यदि (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME) - 1);
    *Info = - info;
    प्रत्यागम;
  }

  तर्काः.alpha = रिक्त;
  तर्काः.beta  = रिक्त;

  *Info = 0;

  यदि (तर्काः.m == 0) प्रत्यागम;

  यदि (diag) {
    यदि (AMIN_K(तर्काः.m, तर्काः.a, तर्काः.lda + 1) == ZERO) {
      *Info = IAMIN_K(तर्काः.m, तर्काः.a, तर्काः.lda + 1);
      प्रत्यागम;
    }
  }


  IDEBUG_START;

  FUNCTION_PROFILE_START();

#यदि_अपरिभाषितम् PPC440
  बफर = (FLOAT *)blas_स्मृतिः_alloc(1);

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_विस्थापनम्_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_विस्थापनम्_B);
#यदि_अन्तः

#यदि_परिभाषितम् SMP
  तर्काः.common = रिक्त;
  यदि (तर्काः.m < 64) {
    तर्काः.nthreads = 1;
  } अन्यथा {
    तर्काः.nthreads = अङ्कः_cpu_avail(4);
  }

  यदि (तर्काः.nthreads == 1) {
#यदि_अन्तः

      (trtrs_single[(uplo << 2) | (trans << 1) | diag])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);

#यदि_परिभाषितम् SMP
  } अन्यथा {
    (trtrs_parallel[(uplo << 2) | (trans << 1) | diag])(&तर्काः, रिक्त, रिक्त, sa, sb, 0);
  }
#यदि_अन्तः

#यदि_अपरिभाषितम् PPC440
  blas_स्मृतिः_free(बफर);
#यदि_अन्तः

  FUNCTION_PROFILE_अन्तः(COMPSIZE * COMPSIZE, तर्काः.m * तर्काः.n, 2 * तर्काः.m * तर्काः.m * तर्काः.n);

  IDEBUG_अन्तः;

  प्रत्यागम;

}
