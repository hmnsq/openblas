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


#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(blasint *N, શૂન્ય *VALPHA, FLOAT *x, blasint *INCX, શૂન્ય *VBETA, FLOAT *y, blasint *INCY)
{

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;
  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;

#નહીંતર

શૂન્ય CNAME(blasint n, શૂન્ય *VALPHA, શૂન્ય *vx, blasint incx, શૂન્ય *VBETA, શૂન્ય *vy, blasint incy)
{
  FLOAT* ALPHA = (FLOAT*) VALPHA;
  FLOAT* BETA = (FLOAT*) VBETA;
  FLOAT* x = (FLOAT*) vx;
  FLOAT* y = (FLOAT*) vy;

#endif

  FLOAT alpha_r = *(ALPHA + 0);
  FLOAT alpha_i = *(ALPHA + 1);
  FLOAT beta_r  = *(BETA  + 0);
  FLOAT beta_i  = *(BETA  + 1);

  જો (n <= 0) પાછા;

  FUNCTION_PROFILE_START();

  જો (incx < 0) x -= (n - 1) * incx * 2;
  જો (incy < 0) y -= (n - 1) * incy * 2;

  AXPBY_K (n, alpha_r, alpha_i, x, incx, beta_r, beta_i, y, incy);

  FUNCTION_PROFILE_END(4, 2 * n, 2 * n);

  પાછા;

}
