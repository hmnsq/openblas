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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(THUNDERX2T99) || defined(VULCAN) || defined(ARMV8) || defined(THUNDERX3T110)
// Multithreaded swap gives performance benefits in ThunderX2T99
#ಇಲ್ಲದಿದ್ದರೆ
// Disable multi-threading as it does not show any performance
// benefits. Keep the multi-threading code for the record.
#ವ್ಯಾಖ್ಯಾನ_ರದ್ದು SMP
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY){

  blasint n    = *N;
  blasint incx = *INCX;
  blasint incy = *INCY;

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx, ಶೂನ್ಯ *vy, blasint incy){
FLOAT *x = (FLOAT*)vx;
FLOAT *y = (FLOAT*)vy;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ;
  FLOAT dummyalpha[2] = {ZERO, ZERO};
  ಪೂರ್ಣಾಂಕ nthreads;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
  PRINT_DEBUG_NAME;
#ಇಲ್ಲದಿದ್ದರೆ
  PRINT_DEBUG_CNAME;
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (n - 1) * incx * 2;
  ಆದರೆ (incy < 0) y -= (n - 1) * incy * 2;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  //disable multi-thread when incx==0 or incy==0
  //In that case, the threads would be dependent.
  ಆದರೆ (incx == 0 || incy == 0 || n < 1048576 * GEMM_MULTITHREAD_ಮಿತಿ / ಗಾತ್ರ(FLOAT))
	  nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ
	  nthreads = ಅಂಕೆ_cpu_avail(1);

  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

  ಅದಲುಬದಲು_K(n, 0, 0, ZERO, ZERO, x, incx, y, incy, ಖಾಲಿ, 0);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
    ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
    ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
    ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ

    blas_level1_thread(ಬಹುಲಕ, n, 0, 0, dummyalpha,
		       x, incx, y, incy, ಖಾಲಿ, 0, (ಶೂನ್ಯ *)ಅದಲುಬದಲು_K, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(2, 2 * n, 0);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;

}
