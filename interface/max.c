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

#ವ್ಯಾಖ್ಯಾನ_ರದ್ದು MAX_K

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ USE_ABS

#ಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ABS fabs
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ABS fabsf
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_MIN

/* ABS & MAX */
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	QAMAX_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	DAMAX_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	SAMAX_K
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	XAMAX_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	ZAMAX_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	CAMAX_K
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಇಲ್ಲದಿದ್ದರೆ

/* ABS & MIN */
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	QAMIN_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	DAMIN_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	SAMIN_K
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	XAMIN_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	ZAMIN_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	CAMIN_K
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ

#ಇಲ್ಲದಿದ್ದರೆ

#ವ್ಯಾಖ್ಯಾನಿಸು ABS

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ USE_MIN

/* MAX */
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	QMAX_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	DMAX_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	SMAX_K
#ಆದರೆ_ಅಂತ್ಯ

#ಇಲ್ಲದಿದ್ದರೆ

/* MIN */
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	QMIN_K
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	DMIN_K
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು MAX_K	SMIN_K
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ 0;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ಆದರೆ (incx == 0) ಹಿಂತಿರುಗಿ (ABS(*x));
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (incx == 0) ಹಿಂತಿರುಗಿ (ABS(*x) + ABS(*(x+1)));
#ಆದರೆ_ಅಂತ್ಯ

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE, n, 0);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ ret;
}

#ಇಲ್ಲದಿದ್ದರೆ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ COMPLEX
FLOAT CNAME(blasint n, ಶೂನ್ಯ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ಇಲ್ಲದಿದ್ದರೆ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#ಆದರೆ_ಅಂತ್ಯ
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ಆದರೆ (n <= 0) ಹಿಂತಿರುಗಿ 0;

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
  ಆದರೆ (incx == 0) ಹಿಂತಿರುಗಿ (ABS(*x));
#ಇಲ್ಲದಿದ್ದರೆ
  ಆದರೆ (incx == 0) ಹಿಂತಿರುಗಿ (ABS(*x) + ABS(*(x+1)));
#ಆದರೆ_ಅಂತ್ಯ

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_ಕೊನೆ(COMPSIZE, n, 0);

  IDEBUG_ಕೊನೆ;

  ಹಿಂತಿರುಗಿ ret;
}

#ಆದರೆ_ಅಂತ್ಯ
