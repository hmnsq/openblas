/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project.                              */
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
  FLOAT alpha = *ALPHA;

#లేకపోతే

శూన్య CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  పూర్ణాంకం బహుళకం, nthreads;
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే CBLAS
  PRINT_DEBUG_NAME;
#లేకపోతే
  PRINT_DEBUG_CNAME;
#అయితే_ముగింపు

  ఒకవేళ (incx <= 0 || n <= 0) తిరిగి;

#నిర్వచితమైతే BGEMM
  దశాంశం alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#లేకపోతే
  FLOAT alpha_float = alpha;
#అయితే_ముగింపు

  ఒకవేళ (alpha_float == ONE) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();


#నిర్వచితమైతే SMP
  ఒకవేళ (n <= 1048576 )
	nthreads = 1;
  లేకపోతే
	nthreads = అంకె_cpu_avail(1);


  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

  SCAL_K(n, 0, 0, alpha, x, incx, ఖాళీ, 0, ఖాళీ, 1);

#నిర్వచితమైతే SMP
  } లేకపోతే {

#నిర్వచితమైతే DOUBLE
    బహుళకం  =  BLAS_DOUBLE | BLAS_REAL;
#లేకపోతే
    బహుళకం  =  BLAS_SINGLE | BLAS_REAL;
#అయితే_ముగింపు

    blas_level1_thread(బహుళకం, n, 0, 0,
#నిర్వచితం_కాకపోతే CBLAS
		       ALPHA,
#లేకపోతే
		       &alpha,
#అయితే_ముగింపు
		       x, incx, ఖాళీ, 0, ఖాళీ, 1,  (పూర్ణాంకం (*)(శూన్య))SCAL_K, nthreads);

  }
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(1, n, n);

  IDEBUG_చివర;

  తిరిగి;

}
