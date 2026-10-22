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

/******************************************************************
 2014/06/07 Saar
******************************************************************/


#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

வெறுமை NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX, FLOAT *BETA, FLOAT *y, blasint *INCY)
{

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT alpha   = *ALPHA;
  FLOAT beta    = *BETA;

#இல்லையெனில்

வெறுமை CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx, FLOAT beta, FLOAT *y, blasint incy)
{

#endif

  எனில் (n <= 0) திரும்பு;

  FUNCTION_PROFILE_START();

  எனில் (incx < 0) x -= (n - 1) * incx;
  எனில் (incy < 0) y -= (n - 1) * incy;

  AXPBY_K(n, alpha, x, incx, beta, y, incy);

  FUNCTION_PROFILE_END(1, 2 * n, 2 * n);

  திரும்பு;

}


