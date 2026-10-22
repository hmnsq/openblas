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

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RETURN_BY_STRUCT
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MYTYPE myxcomplex_t
#ಅಥವಾಆದರೆ defined DOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MYTYPE myzcomplex_t
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MYTYPE myccomplex_t
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RETURN_BY_STRUCT
MYTYPE         NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ಅಥವಾಆದರೆ defined RETURN_BY_ರಾಶಿ
ಶೂನ್ಯ           NAME(OPENBLAS_COMPLEX_FLOAT *ಫಲಿತಾಂಶ, blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ಇಲ್ಲದಿದ್ದರೆ
OPENBLAS_COMPLEX_FLOAT NAME(                        blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY) {
#ಆದರೆ_ಅಂತ್ಯ

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ RETURN_BY_ರಾಶಿ
  OPENBLAS_COMPLEX_FLOAT ret;
#ಆದರೆ_ಅಂತ್ಯ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RETURN_BY_STRUCT
  MYTYPE  myret;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ RETURN_BY_STRUCT
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_NAME;

  ಆದರೆ (n <= 0) {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RETURN_BY_STRUCT
    myret.r = 0.;
    myret.i = 0.;
    ಹಿಂತಿರುಗಿ myret;
#ಅಥವಾಆದರೆ defined RETURN_BY_ರಾಶಿ
    *ಫಲಿತಾಂಶ = zero;
    ಹಿಂತಿರುಗಿ;
#ಇಲ್ಲದಿದ್ದರೆ
    ಹಿಂತಿರುಗಿ zero;
#ಆದರೆ_ಅಂತ್ಯ
  }

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ಆದರೆ (incx < 0) x -= (n - 1) * incx * 2;
  ಆದರೆ (incy < 0) y -= (n - 1) * incy * 2;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ RETURN_BY_STRUCT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಇಲ್ಲದಿದ್ದರೆ
  ret = DOTC_K(n, x, incx, y, incy);
#ಆದರೆ_ಅಂತ್ಯ

  myret.r = CREAL ret;
  myret.i = CIMAG ret;

  FUNCTION_PROFILE_ಕೊನೆ(4, 2 * n, 2 * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ myret;

#ಅಥವಾಆದರೆ defined RETURN_BY_ರಾಶಿ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CONJ
  *ಫಲಿತಾಂಶ = DOTU_K(n, x, incx, y, incy);
#ಇಲ್ಲದಿದ್ದರೆ
  *ಫಲಿತಾಂಶ = DOTC_K(n, x, incx, y, incy);
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(4, 2 * n, 2 * n);

  IDEBUG_ಕೊನೆ;

#ಇಲ್ಲದಿದ್ದರೆ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಇಲ್ಲದಿದ್ದರೆ
  ret = DOTC_K(n, x, incx, y, incy);
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(4, 2 * n, 2 * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ ret;

#ಆದರೆ_ಅಂತ್ಯ

}

#ಇಲ್ಲದಿದ್ದರೆ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FORCE_USE_ರಾಶಿ
ಶೂನ್ಯ           CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx, ಶೂನ್ಯ *vy, blasint incy, ಶೂನ್ಯ* vresult){
OPENBLAS_COMPLEX_FLOAT *ಫಲಿತಾಂಶ= (OPENBLAS_COMPLEX_FLOAT*)vresult;
#ಇಲ್ಲದಿದ್ದರೆ
OPENBLAS_COMPLEX_FLOAT CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx, ಶೂನ್ಯ *vy, blasint incy){

  OPENBLAS_COMPLEX_FLOAT ret;
  OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
#ಆದರೆ_ಅಂತ್ಯ
  FLOAT *x = (FLOAT*) vx;
  FLOAT *y = (FLOAT*) vy;

  PRINT_DEBUG_CNAME;

  ಆದರೆ (n <= 0) {
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FORCE_USE_ರಾಶಿ
    OPENBLAS_COMPLEX_FLOAT zero=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0, 0.0);
    *ಫಲಿತಾಂಶ = zero;
//	CREAL(*result) = 0.0;
//	CIMAG(*result) = 0.0;
    ಹಿಂತಿರುಗಿ;
#ಇಲ್ಲದಿದ್ದರೆ
    ಹಿಂತಿರುಗಿ zero;
#ಆದರೆ_ಅಂತ್ಯ
  }

  ಆದರೆ (incx < 0) x -= (n - 1) * incx * 2;
  ಆದರೆ (incy < 0) y -= (n - 1) * incy * 2;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FORCE_USE_ರಾಶಿ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CONJ
  *ಫಲಿತಾಂಶ = DOTU_K(n, x, incx, y, incy);
#ಇಲ್ಲದಿದ್ದರೆ
  *ಫಲಿತಾಂಶ = DOTC_K(n, x, incx, y, incy);
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(4, 2 * n, 2 * n);

  IDEBUG_ಕೊನೆ;

#ಇಲ್ಲದಿದ್ದರೆ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CONJ
  ret = DOTU_K(n, x, incx, y, incy);
#ಇಲ್ಲದಿದ್ದರೆ
  ret = DOTC_K(n, x, incx, y, incy);
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(4, 2 * n, 2 * n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ ret;

#ಆದರೆ_ಅಂತ್ಯ

}

#ಆದರೆ_ಅಂತ್ಯ
