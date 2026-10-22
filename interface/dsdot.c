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
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

ದ್ವಿಗುಣ NAME(blasint *N, ಭಿನ್ನ *x, blasint *INCX, ಭಿನ್ನ *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  ದ್ವಿಗುಣ ret = 0.0;

  PRINT_DEBUG_NAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (n - 1) * incx;
  ಯದಿ (incy < 0) y -= (n - 1) * incy;

  ret=DSDOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;

}

#ಅನ್ಯಥಾ

ದ್ವಿಗುಣ CNAME(blasint n, ಭಿನ್ನ *x, blasint incx, ಭಿನ್ನ *y, blasint incy){

  ದ್ವಿಗುಣ ret = 0.0;

  PRINT_DEBUG_CNAME;

  ಯದಿ (n <= 0) ಪ್ರತ್ಯಾಗಮ 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ಯದಿ (incx < 0) x -= (n - 1) * incx;
  ಯದಿ (incy < 0) y -= (n - 1) * incy;

  ret=DSDOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ ret;

}

#endif
