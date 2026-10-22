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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGESV"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGESV"
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGESV"
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGESV"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGESV"
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGESV"
#అయితే_ముగింపు
#అయితే_ముగింపు

OPENBLAS_EXPORT
శూన్య NAME(blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA, blasint *ipiv,
         FLOAT *b, blasint *ldB, blasint *Info){

  blas_arg_t వాదనలు;

  blasint info;
  FLOAT *బఫరు;
#నిర్వచితమైతే PPC440
  బాహ్య
#అయితే_ముగింపు
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  వాదనలు.m    = *N;
  వాదనలు.n    = *NRHS;
  వాదనలు.a    = (శూన్య *)a;
  వాదనలు.lda  = *ldA;
  వాదనలు.b    = (శూన్య *)b;
  వాదనలు.ldb  = *ldB;
  వాదనలు.c    = (శూన్య *)ipiv;

  info  = 0;
  ఒకవేళ (వాదనలు.ldb < MAX(1,వాదనలు.m)) info = 7;
  ఒకవేళ (వాదనలు.lda < MAX(1,వాదనలు.m)) info = 4;
  ఒకవేళ (వాదనలు.n   < 0)             info = 2;
  ఒకవేళ (వాదనలు.m   < 0)             info = 1;

  ఒకవేళ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME) - 1);
    *Info = - info;
    తిరిగి;
  }

  వాదనలు.alpha = ఖాళీ;
  వాదనలు.beta  = ఖాళీ;

  *Info = 0;

  ఒకవేళ (వాదనలు.m == 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#నిర్వచితం_కాకపోతే PPC440
  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);
  
  sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  వాదనలు.common = ఖాళీ;

#ఒకవేళ defined(_WIN64) && defined(_M_ARM64)
  #నిర్వచితమైతే COMPLEX
    ఒకవేళ (వాదనలు.m * వాదనలు.m <= 300) 
  #లేకపోతే
    ఒకవేళ (వాదనలు.m * వాదనలు.m <= 500) 
  #అయితే_ముగింపు
      వాదనలు.nthreads = 1;
  లేకపోతే ఒకవేళ (వాదనలు.m * వాదనలు.m <= 1000)
    వాదనలు.nthreads = 4;
  లేకపోతే
    వాదనలు.nthreads = అంకె_cpu_avail(4);
#లేకపోతే
  #నిర్వచితం_కాకపోతే DOUBLE
    ఒకవేళ (వాదనలు.m * వాదనలు.m < 40000)
  #లేకపోతే
    ఒకవేళ (వాదనలు.m * వాదనలు.m < 10000)
  #అయితే_ముగింపు
      వాదనలు.nthreads = 1;
    లేకపోతే
      వాదనలు.nthreads = అంకె_cpu_avail(4);
#అయితే_ముగింపు

  ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    వాదనలు.n    = *N;
    info = GETRF_SINGLE(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

    ఒకవేళ (info == 0 && *NRHS >0){
      వాదనలు.n    = *NRHS;
      GETRS_N_SINGLE(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);
    }

#నిర్వచితమైతే SMP
  } లేకపోతే {

    వాదనలు.n    = *N;
    info = GETRF_PARALLEL(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

    ఒకవేళ (info == 0 && *NRHS > 0){
      వాదనలు.n    = *NRHS;
      GETRS_N_PARALLEL(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);
    }
  }
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే PPC440
  blas_జ్ఞాపకం_free(బఫరు);
#అయితే_ముగింపు

  *Info = info;

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE, *N * *N, 2. / 3. * *N * *N * *N + *N * *N);

  IDEBUG_చివర;

  తిరిగి;
}
