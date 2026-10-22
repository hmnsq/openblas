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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

#undef MAX_K

#ifdef USE_ABS

#ifndef USE_MIN

/* ABS & MAX */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	IQAMAX_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IDAMAX_K
#अण्णहा
#define MAX_K	ISAMAX_K
#endif
#अण्णहा
#ifdef XDOUBLE
#define MAX_K	IXAMAX_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IZAMAX_K
#अण्णहा
#define MAX_K	ICAMAX_K
#endif
#endif

#अण्णहा

/* ABS & MIN */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	IQAMIN_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IDAMIN_K
#अण्णहा
#define MAX_K	ISAMIN_K
#endif
#अण्णहा
#ifdef XDOUBLE
#define MAX_K	IXAMIN_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IZAMIN_K
#अण्णहा
#define MAX_K	ICAMIN_K
#endif
#endif

#endif

#अण्णहा

#ifndef USE_MIN

/* MAX */
#ifdef XDOUBLE
#define MAX_K	IQMAX_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IDMAX_K
#अण्णहा
#define MAX_K	ISMAX_K
#endif

#अण्णहा

/* MIN */
#ifdef XDOUBLE
#define MAX_K	IQMIN_K
#नोचेद्यदि defined(DOUBLE)
#define MAX_K	IDMIN_K
#अण्णहा
#define MAX_K	ISMIN_K
#endif

#endif

#endif

OPENBLAS_EXPORT

#ifndef CBLAS

blasint NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  blasint ret;

  PRINT_DEBUG_NAME;

  जइ (n <= 0) प्रत्यागम 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (blasint)MAX_K(n, x, incx);

  जइ(ret > n) ret=n;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  प्रत्यागम ret;
}

#अण्णहा
#ifdef COMPLEX
CBLAS_INDEX CNAME(blasint n, शून्य *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अण्णहा
CBLAS_INDEX CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  CBLAS_INDEX ret;

  PRINT_DEBUG_CNAME;

  जइ (n <= 0) प्रत्यागम 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  जइ (ret > n) ret=n;

  जइ (ret) ret --;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  प्रत्यागम ret;
}

#endif
