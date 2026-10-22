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
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(blasint *N, FLOAT *ALPHA, FLOAT *x, blasint *INCX){

  blasint n    = *N;
  blasint incx = *INCX;
  FLOAT alpha = *ALPHA;

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(blasint n, FLOAT alpha, FLOAT *x, blasint incx){

#endif

#ifdef SMP
  ಪೂರ್ಣಾಙ್ಕ mode, nthreads;
#endif

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ಅನ್ಯಥಾ
  PRINT_DEBUG_CNAME;
#endif

  ಯದಿ (incx <= 0 || n <= 0) ಪ್ರತ್ಯಾಗಮ;

#ifdef BGEMM
  ಭಿನ್ನ alpha_float;
  SBF16TOS_K(1, &alpha, 1, &alpha_float, 1);
#ಅನ್ಯಥಾ
  FLOAT alpha_float = alpha;
#endif

  ಯದಿ (alpha_float == ONE) ಪ್ರತ್ಯಾಗಮ;

  IDEBUG_START;

  FUNCTION_PROFILE_START();


#ifdef SMP
  ಯದಿ (n <= 1048576 )
	nthreads = 1;
  ಅನ್ಯಥಾ
	nthreads = ಅಂಕೆ_cpu_avail(1);


  ಯದಿ (nthreads == 1) {
#endif

  SCAL_K(n, 0, 0, alpha, x, incx, ರಿಕ್ತ, 0, ರಿಕ್ತ, 1);

#ifdef SMP
  } ಅನ್ಯಥಾ {

#ifdef DOUBLE
    mode  =  BLAS_DOUBLE | BLAS_REAL;
#ಅನ್ಯಥಾ
    mode  =  BLAS_SINGLE | BLAS_REAL;
#endif

    blas_level1_thread(mode, n, 0, 0,
#ifndef CBLAS
		       ALPHA,
#ಅನ್ಯಥಾ
		       &alpha,
#endif
		       x, incx, ರಿಕ್ತ, 0, ರಿಕ್ತ, 1,  (ಪೂರ್ಣಾಙ್ಕ (*)(ಶೂನ್ಯ))SCAL_K, nthreads);

  }
#endif

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;

}
