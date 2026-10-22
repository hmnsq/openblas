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

#ضمن <stdio.h>
#ضمن "common.h"
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#وصف_رد MAX_K

#جيڪڏهن_وصف_ٿيل USE_ABS

#إذا defined(DOUBLE)
#وصف_ڪريو ABS fabs
#وإلا
#وصف_ڪريو ABS fabsf
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل USE_MIN

/* ABS & MAX */
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	QAMAX_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	DAMAX_K
#وإلا
#وصف_ڪريو MAX_K	SAMAX_K
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	XAMAX_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	ZAMAX_K
#وإلا
#وصف_ڪريو MAX_K	CAMAX_K
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#وإلا

/* ABS & MIN */
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	QAMIN_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	DAMIN_K
#وإلا
#وصف_ڪريو MAX_K	SAMIN_K
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	XAMIN_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	ZAMIN_K
#وإلا
#وصف_ڪريو MAX_K	CAMIN_K
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_پورو

#وإلا

#وصف_ڪريو ABS

#جيڪڏهن_وصف_نه_ٿيل USE_MIN

/* MAX */
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	QMAX_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	DMAX_K
#وإلا
#وصف_ڪريو MAX_K	SMAX_K
#جيڪڏهن_پورو

#وإلا

/* MIN */
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو MAX_K	QMIN_K
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو MAX_K	DMIN_K
#وإلا
#وصف_ڪريو MAX_K	SMIN_K
#جيڪڏهن_پورو

#جيڪڏهن_پورو

#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع 0;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#جيڪڏهن_پورو

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE, n, 0);

  IDEBUG_پڇاڙي;

  أرجع ret;
}

#وإلا

#جيڪڏهن_وصف_ٿيل COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#جيڪڏهن_پورو
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0;

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  إذا (incx == 0) أرجع (ABS(*x));
#وإلا
  إذا (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#جيڪڏهن_پورو

  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE, n, 0);

  IDEBUG_پڇاڙي;

  أرجع ret;
}

#جيڪڏهن_پورو
