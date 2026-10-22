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

#ఒకవేళ defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#లేకపోతే
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#నిర్వచనం_రద్దు SMP
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#లేకపోతే

శూన్య CNAME(blasint n, FLOAT *x, blasint incx, FLOAT *y, blasint incy){

#అయితే_ముగింపు

#నిర్వచితమైతే SMP
  పూర్ణాంకం బహుళకం, nthreads;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే CBLAS
  PRINT_DEBUG_NAME;
#లేకపోతే
  PRINT_DEBUG_CNAME;
#అయితే_ముగింపు

  ఒకవేళ (n <= 0) తిరిగి;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (n - 1) * incx;
  ఒకవేళ (incy < 0) y -= (n - 1) * incy;

#నిర్వచితమైతే SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  ఒకవేళ (incx == 0 || incy == 0 || n < 2097152 * GEMM_MULTITHREAD_పరిమితి / పరిమాణం(FLOAT))
    nthreads = 1;
  లేకపోతే
    nthreads = అంకె_cpu_avail(1);

  ఒకవేళ (nthreads == 1) {
#అయితే_ముగింపు

    మార్పిడి_K(n, 0, 0, ZERO, x, incx, y, incy, ఖాళీ, 0);

#నిర్వచితమైతే SMP
  } లేకపోతే {

#నిర్వచితమైతే XDOUBLE
    బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
    బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
    బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు

    blas_level1_thread(బహుళకం, n, 0, 0, dummyalpha,
		       x, incx, y, incy, ఖాళీ, 0, (శూన్య *)మార్పిడి_K, nthreads);
  }

#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(1, 2 * n, 0);

  IDEBUG_చివర;

  తిరిగి;

}
