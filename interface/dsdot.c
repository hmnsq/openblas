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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাৱ "functable.h"
#endif

OPENBLAS_EXPORT
#ifndef CBLAS

দ্ৱিগুণ NAME(blasint *N, ভিন্ন *x, blasint *INCX, ভিন্ন *y, blasint *INCY){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  দ্ৱিগুণ ret = 0.0;

  PRINT_DEBUG_NAME;

  যদি (n <= 0) ঘূৰাই 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (n - 1) * incx;
  যদি (incy < 0) y -= (n - 1) * incy;

  ret=DSDOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  ঘূৰাই ret;

}

#অন্যথা

দ্ৱিগুণ CNAME(blasint n, ভিন্ন *x, blasint incx, ভিন্ন *y, blasint incy){

  দ্ৱিগুণ ret = 0.0;

  PRINT_DEBUG_CNAME;

  যদি (n <= 0) ঘূৰাই 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  যদি (incx < 0) x -= (n - 1) * incx;
  যদি (incy < 0) y -= (n - 1) * incy;

  ret=DSDOT_K(n, x, incx, y, incy);

  FUNCTION_PROFILE_END(1, n, n);

  IDEBUG_END;

  ঘূৰাই ret;

}

#endif
