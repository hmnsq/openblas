/***************************************************************************
Copyright (c) 2013, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/**************************************************************************************
* 2014/05/02 Saar
* 	fixed two bugs as reported by Brendan Tracey
*	Test with lapack-3.5.0	: OK
*
**************************************************************************************/


#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#నిర్వచించు  GAM     4096.e0
#నిర్వచించు  GAMSQ   16777216.e0
#నిర్వచించు  RGAMSQ  5.9604645e-8

#నిర్వచించు  TWO 2.e0

#నిర్వచితమైతే DOUBLE
#నిర్వచించు ABS(x) fabs(x)
#లేకపోతే
#నిర్వచించు ABS(x) fabsf(x)
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(FLOAT *dd1, FLOAT *dd2, FLOAT *dx1, FLOAT *DY1, FLOAT *dparam){

  FLOAT	dy1 = *DY1;

#లేకపోతే

శూన్య CNAME(FLOAT *dd1, FLOAT *dd2, FLOAT *dx1, FLOAT dy1, FLOAT *dparam){

#అయితే_ముగింపు

	FLOAT du, dp1, dp2, dq2, dq1, dh11=ZERO, dh21=ZERO, dh12=ZERO, dh22=ZERO, dflag=-ONE, dtemp;

	ఒకవేళ (*dd2 == ZERO || dy1 == ZERO)
	{
		dflag = -TWO;
		dparam[0] = dflag;
		తిరిగి;
	}
		
	ఒకవేళ(*dd1 < ZERO)
	{
		dflag = -ONE;
		dh11  = ZERO;
		dh12  = ZERO;
		dh21  = ZERO;
		dh22  = ZERO;

		*dd1  = ZERO;
		*dd2  = ZERO;
		*dx1  = ZERO;
	}
	లేకపోతే ఒకవేళ ((*dd1 == ZERO || *dx1 == ZERO) && *dd2 > ZERO)
	{
		dflag = ONE;
		dh12 = 1;
		dh21 = -1;
		*dx1 = dy1;
		dtemp = *dd1;
		*dd1 = *dd2;
		*dd2 = dtemp;
	} 
	లేకపోతే
	{
		dp2 = *dd2 * dy1;
		dp1 = *dd1 * *dx1;
		dq2 =  dp2 * dy1;
		dq1 =  dp1 * *dx1;
		ఒకవేళ(ABS(dq1) > ABS(dq2))
		{
			dh11  =  ONE;
			dh22  =  ONE;
			dh21 = -  dy1 / *dx1;
			dh12 =    dp2 /  dp1;

			du   = ONE - dh12 * dh21;
			ఒకవేళ(du <= ZERO)
			{
				/* This can only happen through rounding: the
				   branch was chosen on ABS(dq1) > ABS(dq2),
				   which in exact arithmetic makes du positive.
				   Matches the reference BLAS srotmg.f and the
				   gonum implementation. */
				dflag = -ONE;

				dh11  = ZERO;
				dh12  = ZERO;
				dh21  = ZERO;
				dh22  = ZERO;

				*dd1  = ZERO;
				*dd2  = ZERO;
				*dx1  = ZERO;
			}
			లేకపోతే
			{
			dflag = ZERO;
			*dd1  = *dd1 / du;
			*dd2  = *dd2 / du;
			*dx1  = *dx1 * du;
			}
		}
		లేకపోతే
		{
			ఒకవేళ(dq2 < ZERO)
			{
				dflag = -ONE;

				dh11  = ZERO;
				dh12  = ZERO;
				dh21  = ZERO;
				dh22  = ZERO;

				*dd1  = ZERO;
				*dd2  = ZERO;
				*dx1  = ZERO;
			}
			లేకపోతే
			{
				dflag =  ONE;
				dh21  = -ONE;
				dh12  =  ONE;

				dh11  =  dp1 /  dp2;
				dh22  = *dx1 /  dy1;
				du    =  ONE + dh11 * dh22;
				dtemp = *dd2 / du;

				*dd2  = *dd1 / du;
				*dd1  = dtemp;
				*dx1  = dy1 * du;
			}
		}


		అయితే ( *dd1 <= RGAMSQ && *dd1 != ZERO)
		{
			dflag = -ONE;
			*dd1  = *dd1 * (GAM * GAM);
			*dx1  = *dx1 / GAM;
			dh11  = dh11 / GAM;
			dh12  = dh12 / GAM;
		}
		అయితే (ABS(*dd1) > GAMSQ) {
			dflag = -ONE;
			*dd1  = *dd1 / (GAM * GAM);
			*dx1  = *dx1 * GAM;
			dh11  = dh11 * GAM;
			dh12  = dh12 * GAM;
		}

		అయితే (ABS(*dd2) <= RGAMSQ && *dd2 != ZERO) {
			dflag = -ONE;
			*dd2  = *dd2 * (GAM * GAM);
			dh21  = dh21 / GAM;
			dh22  = dh22 / GAM;
		}
		అయితే (ABS(*dd2) > GAMSQ) {
			dflag = -ONE;
			*dd2  = *dd2 / (GAM * GAM);
			dh21  = dh21 * GAM;
			dh22  = dh22 * GAM;
		}

	}

	ఒకవేళ(dflag < ZERO)
	{
		dparam[1] = dh11;
		dparam[2] = dh21;
		dparam[3] = dh12;
		dparam[4] = dh22;
	}
	లేకపోతే
	{
		ఒకవేళ(dflag == ZERO)
		{
			dparam[2] = dh21;
			dparam[3] = dh12;
		}
		లేకపోతే
		{
			dparam[1] = dh11;
			dparam[4] = dh22;
		}
	}


	dparam[0] = dflag;
	తిరిగి;
}


