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

#यदि defined(DOUBLE)
#define ABS भ_असल
#अन्यथा
#define ABS fabsf
#endif

#ifndef USE_MIN

/* ABS & MAX */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMAX_K
#अथवायदि defined(DOUBLE)
#define MAX_K	DAMAX_K
#अन्यथा
#define MAX_K	SAMAX_K
#endif
#अन्यथा
#ifdef XDOUBLE
#define MAX_K	XAMAX_K
#अथवायदि defined(DOUBLE)
#define MAX_K	ZAMAX_K
#अन्यथा
#define MAX_K	CAMAX_K
#endif
#endif

#अन्यथा

/* ABS & MIN */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMIN_K
#अथवायदि defined(DOUBLE)
#define MAX_K	DAMIN_K
#अन्यथा
#define MAX_K	SAMIN_K
#endif
#अन्यथा
#ifdef XDOUBLE
#define MAX_K	XAMIN_K
#अथवायदि defined(DOUBLE)
#define MAX_K	ZAMIN_K
#अन्यथा
#define MAX_K	CAMIN_K
#endif
#endif

#endif

#अन्यथा

#define ABS

#ifndef USE_MIN

/* MAX */
#ifdef XDOUBLE
#define MAX_K	QMAX_K
#अथवायदि defined(DOUBLE)
#define MAX_K	DMAX_K
#अन्यथा
#define MAX_K	SMAX_K
#endif

#अन्यथा

/* MIN */
#ifdef XDOUBLE
#define MAX_K	QMIN_K
#अथवायदि defined(DOUBLE)
#define MAX_K	DMIN_K
#अन्यथा
#define MAX_K	SMIN_K
#endif

#endif

#endif

OPENBLAS_EXPORT

#ifndef CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  यदि (n <= 0) फर्काऊ 0;

#ifndef COMPLEX
  यदि (incx == 0) फर्काऊ (ABS(*x));
#अन्यथा
  यदि (incx == 0) फर्काऊ (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  फर्काऊ ret;
}

#अन्यथा

#ifdef COMPLEX
FLOAT CNAME(blasint n, शून्य *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  यदि (n <= 0) फर्काऊ 0;

#ifndef COMPLEX
  यदि (incx == 0) फर्काऊ (ABS(*x));
#अन्यथा
  यदि (incx == 0) फर्काऊ (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  फर्काऊ ret;
}

#endif
