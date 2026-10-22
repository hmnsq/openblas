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

/**********************************************************************
 2014/06/07 Saar

**********************************************************************/


#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

व्योम NAME(blasint *N, व्योम *VALPHA, FLOAT *x, blasint *INCX, व्योम *VBETA, FLOAT *y, blasint *INCY)
{

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;

#अन्यथा

व्योम CNAME(blasint n, व्योम *VALPHA, व्योम *vx, blasint incx, व्योम *VBETA, व्योम *vy, blasint incy)
{
  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

#पूर्ण_अगर

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);
  FLOAT beta_r  = *(BETA  + 0);
  FLOAT beta_i  = *(BETA  + 1);

  अगर (n <= 0) वापस;

  FUNCTION_PROFILE_START();

  अगर (incx < 0) x -= (n - 1) * incx * 2;
  अगर (incy < 0) y -= (n - 1) * incy * 2;

  AXPBY_K (n, alpha_r, alpha_i, x, incx, beta_r, beta_i, y, incy);

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  वापस;

}
