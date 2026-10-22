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

#undef MAX_K

#ifdef USE_ABS

#ifndef USE_MIN

/* ABS & MAX */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	IQAMAX_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IDAMAX_K
#অন্যথা
#define MAX_K	ISAMAX_K
#endif
#অন্যথা
#ifdef XDOUBLE
#define MAX_K	IXAMAX_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IZAMAX_K
#অন্যথা
#define MAX_K	ICAMAX_K
#endif
#endif

#অন্যথা

/* ABS & MIN */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	IQAMIN_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IDAMIN_K
#অন্যথা
#define MAX_K	ISAMIN_K
#endif
#অন্যথা
#ifdef XDOUBLE
#define MAX_K	IXAMIN_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IZAMIN_K
#অন্যথা
#define MAX_K	ICAMIN_K
#endif
#endif

#endif

#অন্যথা

#ifndef USE_MIN

/* MAX */
#ifdef XDOUBLE
#define MAX_K	IQMAX_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IDMAX_K
#অন্যথা
#define MAX_K	ISMAX_K
#endif

#অন্যথা

/* MIN */
#ifdef XDOUBLE
#define MAX_K	IQMIN_K
#নোচেদ্যদি defined(DOUBLE)
#define MAX_K	IDMIN_K
#অন্যথা
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

  যদি (n <= 0) ঘূৰাই 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (blasint)MAX_K(n, x, incx);

  যদি(ret > n) ret=n;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  ঘূৰাই ret;
}

#অন্যথা
#ifdef COMPLEX
CBLAS_INDEX CNAME(blasint n, শূন্য *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#অন্যথা
CBLAS_INDEX CNAME(blasint n, FLOAT *x, blasint incx){
#endif

  CBLAS_INDEX ret;

  PRINT_DEBUG_CNAME;

  যদি (n <= 0) ঘূৰাই 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  যদি (ret > n) ret=n;

  যদি (ret) ret --;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  ঘূৰাই ret;
}

#endif
