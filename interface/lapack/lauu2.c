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
#నిర్వచించు ERROR_NAME "QLAUU2"
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DLAUU2"
#లేకపోతే
#నిర్వచించు ERROR_NAME "SLAUU2"
#అయితే_ముగింపు

స్థిర blasint (*lauu2[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#నిర్వచితమైతే XDOUBLE
  qlauu2_U, qlauu2_L,
#లేదావేళ defined(DOUBLE)
  dlauu2_U, dlauu2_L,
#లేకపోతే
  slauu2_U, slauu2_L,
#అయితే_ముగింపు
  };

OPENBLAS_EXPORT
శూన్య NAME(అక్షరం *UPLO, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t వాదనలు;

  blasint uplo_arg = *UPLO;
  blasint uplo;
  blasint info;
  FLOAT *బఫరు;
#నిర్వచితమైతే PPC440
  బాహ్య
#అయితే_ముగింపు
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  వాదనలు.n    = *N;
  వాదనలు.a    = (శూన్య *)a;
  వాదనలు.lda  = *ldA;

  TOUPPER(uplo_arg);

  uplo = -1;
  ఒకవేళ (uplo_arg == 'U') uplo = 0;
  ఒకవేళ (uplo_arg == 'L') uplo = 1;

  info  = 0;
  ఒకవేళ (వాదనలు.lda < MAX(1,వాదనలు.n)) info = 4;
  ఒకవేళ (వాదనలు.n   < 0)             info = 2;
  ఒకవేళ (uplo     < 0)             info = 1;
  ఒకవేళ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME) - 1);
    *Info = - info;
    తిరిగి;
  }

  *Info = 0;

  ఒకవేళ (వాదనలు.n <= 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#నిర్వచితం_కాకపోతే PPC440
  బఫరు = (FLOAT *)blas_జ్ఞాపకం_alloc(1);

  sa = (FLOAT *)((BLASLONG)బఫరు + GEMM_స్థానభ్రంశం_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * కొలత + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_స్థానభ్రంశం_B);
#అయితే_ముగింపు

  info = (lauu2[uplo])(&వాదనలు, ఖాళీ, ఖాళీ, sa, sb, 0);

  *Info = info;

#నిర్వచితం_కాకపోతే PPC440
  blas_జ్ఞాపకం_free(బఫరు);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(1, .5 * వాదనలు.n * వాదనలు.n,
		       వాదనలు.n * (1./3. + వాదనలు.n * ( 1./2. + వాదనలు.n * 1./6.))
		       +  1./6. * వాదనలు.n * (వాదనలు.n * వాదనలు.n - 1));

  IDEBUG_చివర;

  తిరిగి;
}
