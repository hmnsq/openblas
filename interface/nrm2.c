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

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ఒకవేళ (n <= 0) తిరిగి 0.;

#నిర్వచితం_కాకపోతే COMPLEX
  ఒకవేళ (n == 1)
#నిర్వచితమైతే DOUBLE
    తిరిగి fabs(x[0]);
#లేకపోతే
    తిరిగి fabsf(x[0]);
#అయితే_ముగింపు
#అయితే_ముగింపు

  ఒకవేళ (incx == 0)
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే DOUBLE
  తిరిగి (వర్గమూలం((ద్విగుణ)n)*fabs(x[0]));
#లేకపోతే
  తిరిగి (వర్గమూలం((దశాంశం)n)*fabsf(x[0]));
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే DOUBLE
  {
  ద్విగుణ fr=fabs(x[0]);
  ద్విగుణ fi=fabs(x[1]);
  ద్విగుణ fmin=MIN(fr,fi);
  ద్విగుణ fmax=MAX(fr,fi);
  ఒకవేళ (fmax==0.) తిరిగి(fmax);
  ఒకవేళ (fmax==fmin) తిరిగి(వర్గమూలం((ద్విగుణ)n)*వర్గమూలం(2.)*fmax);
  తిరిగి (వర్గమూలం((ద్విగుణ)n) * fmax * వర్గమూలం (1. + (fmin/fmax)*(fmin/fmax)));
  }
#లేకపోతే
  {
  దశాంశం fr=fabs(x[0]);
  దశాంశం fi=fabs(x[1]);
  దశాంశం fmin=MIN(fr,fi);
  దశాంశం fmax=MAX(fr,fi);
  ఒకవేళ (fmax==0.) తిరిగి(fmax);
  ఒకవేళ (fmax==fmin) తిరిగి(వర్గమూలం((దశాంశం)n)*వర్గమూలం(2.)*fmax);
  తిరిగి (వర్గమూలం((దశాంశం)n) * fmax * వర్గమూలం (1. + (fmin/fmax)*(fmin/fmax)));
  }
#అయితే_ముగింపు
#అయితే_ముగింపు

  ఒకవేళ (incx < 0) 
#నిర్వచితమైతే COMPLEX    
    x -= (n - 1) * incx * 2;
#లేకపోతే
    x -= (n - 1) * incx;
#అయితే_ముగింపు
  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ret = (FLOATRET)NRM2_K(n, x, incx);

  FUNCTION_PROFILE_చివర(COMPSIZE, n, 2 * n);

  IDEBUG_చివర;

  తిరిగి ret;
}

#లేకపోతే

#నిర్వచితమైతే COMPLEX
FLOAT CNAME(blasint n, శూన్య *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#లేకపోతే
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#అయితే_ముగింపు

  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ఒకవేళ (n <= 0) తిరిగి 0.;

#నిర్వచితం_కాకపోతే COMPLEX
  ఒకవేళ (n == 1)
#నిర్వచితమైతే DOUBLE
    తిరిగి fabs(x[0]);
#లేకపోతే
    తిరిగి fabsf(x[0]);
#అయితే_ముగింపు
#అయితే_ముగింపు

  ఒకవేళ (incx == 0)
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే DOUBLE
  తిరిగి (వర్గమూలం((ద్విగుణ)n)*fabs(x[0]));
#లేకపోతే
  తిరిగి (వర్గమూలం((దశాంశం)n)*fabsf(x[0]));
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే DOUBLE
  {
  ద్విగుణ fr=fabs(x[0]);
  ద్విగుణ fi=fabs(x[1]);
  ద్విగుణ fmin=MIN(fr,fi);
  ద్విగుణ fmax=MAX(fr,fi);
  ఒకవేళ (fmax==0.) తిరిగి(fmax);
  ఒకవేళ (fmax==fmin) తిరిగి(వర్గమూలం((ద్విగుణ)n)*వర్గమూలం(2.)*fmax);
  తిరిగి (వర్గమూలం((ద్విగుణ)n) * fmax * వర్గమూలం (1. + (fmin/fmax)*(fmin/fmax)));
  }
#లేకపోతే
  {
  దశాంశం fr=fabs(x[0]);
  దశాంశం fi=fabs(x[1]);
  దశాంశం fmin=MIN(fr,fi);
  దశాంశం fmax=MAX(fr,fi);
  ఒకవేళ (fmax==0.) తిరిగి(fmax);
  ఒకవేళ (fmax==fmin) తిరిగి(వర్గమూలం((దశాంశం)n)*వర్గమూలం(2.)*fmax);
  తిరిగి (వర్గమూలం((దశాంశం)n) * fmax * వర్గమూలం (1. + (fmin/fmax)*(fmin/fmax)));
  }
#అయితే_ముగింపు
#అయితే_ముగింపు

  ఒకవేళ (incx < 0) 
#నిర్వచితమైతే COMPLEX    
    x -= (n - 1) * incx * 2;
#లేకపోతే
    x -= (n - 1) * incx;
#అయితే_ముగింపు
  
  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();

  ret = NRM2_K(n, x, incx);

  FUNCTION_PROFILE_చివర(COMPSIZE, n, 2 * n);

  IDEBUG_చివర;

  తిరిగి ret;
}

#అయితే_ముగింపు
