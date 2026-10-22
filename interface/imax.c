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

#నిర్వచనం_రద్దు MAX_K

#నిర్వచితమైతే USE_ABS

#నిర్వచితం_కాకపోతే USE_MIN

/* ABS & MAX */
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IQAMAX_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IDAMAX_K
#లేకపోతే
#నిర్వచించు MAX_K	ISAMAX_K
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IXAMAX_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IZAMAX_K
#లేకపోతే
#నిర్వచించు MAX_K	ICAMAX_K
#అయితే_ముగింపు
#అయితే_ముగింపు

#లేకపోతే

/* ABS & MIN */
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IQAMIN_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IDAMIN_K
#లేకపోతే
#నిర్వచించు MAX_K	ISAMIN_K
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IXAMIN_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IZAMIN_K
#లేకపోతే
#నిర్వచించు MAX_K	ICAMIN_K
#అయితే_ముగింపు
#అయితే_ముగింపు

#అయితే_ముగింపు

#లేకపోతే

#నిర్వచితం_కాకపోతే USE_MIN

/* MAX */
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IQMAX_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IDMAX_K
#లేకపోతే
#నిర్వచించు MAX_K	ISMAX_K
#అయితే_ముగింపు

#లేకపోతే

/* MIN */
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు MAX_K	IQMIN_K
#లేదావేళ defined(DOUBLE)
#నిర్వచించు MAX_K	IDMIN_K
#లేకపోతే
#నిర్వచించు MAX_K	ISMIN_K
#అయితే_ముగింపు

#అయితే_ముగింపు

#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

blasint NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  blasint ret;

  PRINT_DEBUG_NAME;

  ఒకవేళ (n <= 0) తిరిగి 0;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ret = (blasint)MAX_K(n, x, incx);

  ఒకవేళ(ret > n) ret=n;

  FUNCTION_PROFILE_చివర(COMPSIZE, n, 0);

  IDEBUG_చివర;

  తిరిగి ret;
}

#లేకపోతే
#నిర్వచితమైతే COMPLEX
CBLAS_సూచీ CNAME(blasint n, శూన్య *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#లేకపోతే
CBLAS_సూచీ CNAME(blasint n, FLOAT *x, blasint incx){
#అయితే_ముగింపు

  CBLAS_సూచీ ret;

  PRINT_DEBUG_CNAME;

  ఒకవేళ (n <= 0) తిరిగి 0;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ret = MAX_K(n, x, incx);

  ఒకవేళ (ret > n) ret=n;

  ఒకవేళ (ret) ret --;

  FUNCTION_PROFILE_చివర(COMPSIZE, n, 0);

  IDEBUG_చివర;

  తిరిగి ret;
}

#అయితే_ముగింపు
