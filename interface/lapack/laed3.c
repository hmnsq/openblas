/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"

#ਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DLAED3"
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SLAED3"
#ਜੇ_ਅੰਤ

/* ===================================================================== */
OPENBLAS_EXPORT
ਖਾਲੀ NAME(blasint *k, blasint *n, blasint *n1, FLOAT *d, 
        FLOAT *q, blasint *ldq, FLOAT *rho, FLOAT *dlamda,
        FLOAT *q2, blasint *indx, blasint *ctot, FLOAT *w, 
        FLOAT *s, blasint *Info)
{
  blasint kval, nval, qdim, info;

  qdim = *ldq;
  kval = *k;
  nval = *n;

/*   Test the input parameters. */
  info = 0;
  ਜੇ (kval < 0) {
    info = 1;
  } ਨਹੀਂਤਾਂ ਜੇ (nval < kval) {
    info = 2;
  } ਨਹੀਂਤਾਂ ਜੇ (qdim < nval || qdim < 1) {
    info = 6;
  }
  ਜੇ (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME) - 1);
    *Info = - info;
    ਵਾਪਸ;
  }

/*   Quick return if possible */

  *Info = 0;
  ਜੇ (kval == 0) ਵਾਪਸ;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਪੂਰਨ nthreads = 1;
  ਜੇ (nval >= 64) nthreads = ਅੰਕ_cpu_avail(4);

  ਜੇ (nthreads == 1) {
#ਜੇ_ਅੰਤ
    LAED3_SINGLE(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, Info);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {
    LAED3_PARALLEL(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, Info);
  }
#ਜੇ_ਅੰਤ

  ਵਾਪਸ;
}
