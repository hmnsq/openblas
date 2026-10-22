/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project.                              */
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

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ, nthreads;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
  PRINT_DEBUG_NAME;
#ಇಲ್ಲದಿದ್ದರೆ
  PRINT_DEBUG_CNAME;
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (incx <= 0 || n <= 0) ಹಿಂತಿರುಗಿ;

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ BGEMM
  ದಶಮಾಂಶ alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#ಇಲ್ಲದಿದ್ದರೆ
  FLOAT alpha_float = alpha;
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (alpha_float == ONE) ಹಿಂತಿರುಗಿ;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();


#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  ಆದರೆ (n <= 1048576 )
	nthreads = 1;
  ಇಲ್ಲದಿದ್ದರೆ
	nthreads = ಅಂಕೆ_cpu_avail(1);


  ಆದರೆ (nthreads == 1) {
#ಆದರೆ_ಅಂತ್ಯ

  SCAL_K(n, 0, 0, alpha, x, incx, ಖಾಲಿ, 0, ಖಾಲಿ, 1);

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
  } ಇಲ್ಲದಿದ್ದರೆ {

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DOUBLE
    ಬಹುಲಕ  =  BLAS_DOUBLE | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
    ಬಹುಲಕ  =  BLAS_SINGLE | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ

    blas_level1_thread(ಬಹುಲಕ, n, 0, 0,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
		       ALPHA,
#ಇಲ್ಲದಿದ್ದರೆ
		       &alpha,
#ಆದರೆ_ಅಂತ್ಯ
		       x, incx, ಖಾಲಿ, 0, ಖಾಲಿ, 1,  (ಪೂರ್ಣಾಂಕ (*)(ಶೂನ್ಯ))SCAL_K, nthreads);

  }
#ಆದರೆ_ಅಂತ್ಯ

  FUNCTION_PROFILE_ಕೊನೆ(1, n, n);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ;

}
