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

#నిర్వచితమైతే RETURN_BY_STRUCT
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MYTYPE myxcomplex_t
#లేదావేళ defined DOUBLE
#నిర్వచించు MYTYPE myzcomplex_t
#లేకపోతే
#నిర్వచించు MYTYPE myccomplex_t
#అయితే_ముగింపు
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

#నిర్వచితమైతే RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#లేదావేళ defined RETURN_BY_దొంతర
శూన్య           NAME(OPENBLAS_COMPLEX_FLOAT *ఫలితం, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#లేకపోతే
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#అయితే_ముగింపు

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#నిర్వచితం_కాకపోతే RETURN_BY_దొంతర
  OPENBLAS_COMPLEX_FLOAT ret;
#అయితే_ముగింపు
#నిర్వచితమైతే RETURN_BY_STRUCT
  MYTYPE  myret;
#అయితే_ముగింపు

#నిర్వచితం_కాకపోతే RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#అయితే_ముగింపు

  PRINT_DEBUG_NAME;

  ఒకవేళ (n <= 0) {
#నిర్వచితమైతే RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    తిరిగి myret;
#లేదావేళ defined RETURN_BY_దొంతర
    *ఫలితం = zero;
    తిరిగి;
#లేకపోతే
    తిరిగి zero;
#అయితే_ముగింపు
  }

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (incx < 0) x -= (n - 1) * incx * 2;
  ఒకవేళ (incy < 0) y -= (n - 1) * incy * 2;

#నిర్వచితమైతే RETURN_BY_STRUCT

#నిర్వచితం_కాకపోతే CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#లేకపోతే
  ret = DOTC_K(n, x, incx, y, incy);
#అయితే_ముగింపు

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_చివర(4, 2 * n, 2 * n);

  IDEBUG_చివర;

  తిరిగి myret;

#లేదావేళ defined RETURN_BY_దొంతర

#నిర్వచితం_కాకపోతే CONJ
  *ఫలితం = DOTU_K(n, x, incx, y, incy);
#లేకపోతే
  *ఫలితం = DOTC_K(n, x, incx, y, incy);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(4, 2 * n, 2 * n);

  IDEBUG_చివర;

#లేకపోతే

#నిర్వచితం_కాకపోతే CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#లేకపోతే
  ret = DOTC_K(n, x, incx, y, incy);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(4, 2 * n, 2 * n);

  IDEBUG_చివర;

  తిరిగి ret;

#అయితే_ముగింపు

}

#లేకపోతే

#నిర్వచితమైతే FORCE_USE_దొంతర
శూన్య           CNAME(blasint n, శూన్య *vx, blasint incx, శూన్య *vy, blasint incy, శూన్య* vresult){
OPENBLAS_COMPLEX_FLOAT *ఫలితం= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#లేకపోతే
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, శూన్య *vx, blasint incx, శూన్య *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#అయితే_ముగింపు
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  ఒకవేళ (n <= 0) {
#నిర్వచితమైతే FORCE_USE_దొంతర
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *ఫలితం = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    తిరిగి;
#లేకపోతే
    తిరిగి zero;
#అయితే_ముగింపు
  }

  ఒకవేళ (incx < 0) x -= (n - 1) * incx * 2;
  ఒకవేళ (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

#నిర్వచితమైతే FORCE_USE_దొంతర

#నిర్వచితం_కాకపోతే CONJ
  *ఫలితం = DOTU_K(n, x, incx, y, incy);
#లేకపోతే
  *ఫలితం = DOTC_K(n, x, incx, y, incy);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(4, 2 * n, 2 * n);

  IDEBUG_చివర;

#లేకపోతే

#నిర్వచితం_కాకపోతే CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#లేకపోతే
  ret = DOTC_K(n, x, incx, y, incy);
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(4, 2 * n, 2 * n);

  IDEBUG_చివర;

  తిరిగి ret;

#అయితే_ముగింపు

}

#అయితే_ముగింపు
