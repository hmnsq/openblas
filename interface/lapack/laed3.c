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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"

#यदि defined(DOUBLE)
#define ERROR_NAME "DLAED3"
#अन्यथा
#define ERROR_NAME "SLAED3"
#endif

/* ===================================================================== */
OPENBLAS_EXPORT
शून्य NAME(blasint *k, blasint *n, blasint *n1, FLOAT *d, 
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
  यदि (kval < 0) {
    info = 1;
  } अन्यथा यदि (nval < kval) {
    info = 2;
  } अन्यथा यदि (qdim < nval || qdim < 1) {
    info = 6;
  }
  यदि (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME) - 1);
    *Info = - info;
    प्रत्यागम;
  }

/*   Quick return if possible */

  *Info = 0;
  यदि (kval == 0) प्रत्यागम;

#ifdef SMP
  पूर्णाङ्क nthreads = 1;
  यदि (nval >= 64) nthreads = अंक_cpu_avail(4);

  यदि (nthreads == 1) {
#endif
    LAED3_SINGLE(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, Info);
#ifdef SMP
  } अन्यथा {
    LAED3_PARALLEL(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, Info);
  }
#endif

  प्रत्यागम;
}
