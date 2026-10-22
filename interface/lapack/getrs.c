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

#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGETRS"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGETRS"
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGETRS"
#అయితే_ముగింపు

స్థిర blasint (*getrs_single[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_SINGLE, GETRS_T_SINGLE,
};

#నిర్వచితమైతే SMP
స్థిర blasint (*getrs_parallel[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
  GETRS_N_PARALLEL, GETRS_T_PARALLEL,
};
#అయితే_ముగింపు

OPENBLAS_EXPORT
శూన్య NAME(అక్షరం *TRANS, blasint *N, blasint *NRHS, FLOAT *a, blasint *ldA,
  blasint *ipiv, FLOAT *b, blasint *ldB, blasint *Info){

  అక్షరం trans_arg = *TRANS;

  blas_arg_t వాదనలు;

  blasint info;
  పూర్ణాంకం trans;
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

  info = 0;

  TOUPPER(trans_arg);
  trans = -1;

  ఒకవేళ (trans_arg == 'N') trans = 0;
  ఒకవేళ (trans_arg == 'T') trans = 1;
  ఒకవేళ (trans_arg == 'R') trans = 0;
  ఒకవేళ (trans_arg == 'C') trans = 1;

  ఒకవేళ (వాదనలు.ldb  < MAX(1, వాదనలు.m)) info = 8;
  ఒకవేళ (వాదనలు.lda  < MAX(1, వాదనలు.m)) info = 5;
  ఒకవేళ (వాదనలు.n    < 0) info = 3;
  ఒకవేళ (వాదనలు.m    < 0) info = 2;
  ఒకవేళ (trans     < 0) info = 1;

  ఒకవేళ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME) - 1);
    తిరిగి;
  }

  వాదనలు.alpha = ఖాళీ;
  వాదనలు.beta  = ఖాళీ;

  *Info = info;

  ఒకవేళ (వాదనలు.m == 0 || వాదనలు.n == 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#నిర్వచితం_కాకపోతే PPC440
  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);

  sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  వాదనలు.common = ఖాళీ;
  ఒకవేళ (వాదనలు.m < 64) {
    వాదనలు.nthreads = 1;
  } లేకపోతే {
    వాదనలు.nthreads = అంకె_cpu_avail(4);
  }

  ఒకవేళ (వాదనలు.nthreads == 1) {
#అయితే_ముగింపు

    (getrs_single[trans])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

#నిర్వచితమైతే SMP
  } లేకపోతే {
    (getrs_parallel[trans])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);
  }
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే PPC440
  blas_జ్ఞాపకం_free(బఫరు);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(COMPSIZE * COMPSIZE, వాదనలు.m * వాదనలు.n, 2 * వాదనలు.m * వాదనలు.m * వాదనలు.n);

  IDEBUG_చివర;

  తిరిగి;

}
