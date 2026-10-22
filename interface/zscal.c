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

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;

#నిర్వచితం_కాకపోతే SSCAL
  FLOAT *alpha=ALPHA;
#లేకపోతే
  FLOAT alpha[2] = {ALPHA[0], ZERO};
#అయితే_ముగింపు

#లేకపోతే

#నిర్వచితం_కాకపోతే SSCAL
శూన్య CNAME(blasint n, శూన్య *VALPHA, శూన్య *vx, blasint incx){

  FLOAT *x = (FLOAT*) vx;
  FLOAT *alpha=(FLOAT*)VALPHA;
#లేకపోతే
శూన్య CNAME(blasint n, FLOAT alpha_r, శూన్య *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;

  FLOAT alpha[2] = {alpha_r, ZERO};
#అయితే_ముగింపు
#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  పూర్ణాంకం బహుళకం;
  పూర్ణాంకం nthreads;
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే CBLAS
  PRINT_DEBUG_NAME;
#లేకపోతే
  PRINT_DEBUG_CNAME;
#అయితే_ముగింపు

  ఒకవేళ (incx <= 0 || n <= 0) తిరిగి;

  ఒకవేళ ((alpha[0] == ONE) && (alpha[1] == ZERO)) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#నిర్వచితమైతే SMP
  ఒకవేళ ( n <= 1048576 )
	nthreads = 1;
  లేకపోతే
	nthreads = అంకె_cpu_avail(1);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  SCAL_K(n, 0, 0, alpha[0], alpha[1], x, incx, ఖాళీ, 0, ఖాళీ, 1);

#నిర్వచితమైతే SMP
  } లేకపోతే {
#నిర్వచితమైతే DOUBLE
    బహుళకం  =  BLAS_DOUBLE | BLAS_COMPLEX;
#లేకపోతే
    బహుళకం  =  BLAS_SINGLE | BLAS_COMPLEX;
#అయితే_ముగింపు

    blas_level1_thread(బహుళకం, n, 0, 0,  alpha, x, incx, ఖాళీ, 0, ఖాళీ, 1, (పూర్ణాంకం (*)(శూన్య))SCAL_K, nthreads);

  }
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(4, n, n);

  IDEBUG_చివర;

  తిరిగి;

}
