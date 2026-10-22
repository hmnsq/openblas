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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#undef MAX_K

#ifdef USE_ABS

#𐤀𐤌 defined(DOUBLE)
#define ABS fabs
#وإلا
#define ABS fabsf
#endif

#ifndef USE_MIN

/* ABS & MAX */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMAX_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	DAMAX_K
#وإلا
#define MAX_K	SAMAX_K
#endif
#وإلا
#ifdef XDOUBLE
#define MAX_K	XAMAX_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	ZAMAX_K
#وإلا
#define MAX_K	CAMAX_K
#endif
#endif

#وإلا

/* ABS & MIN */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMIN_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	DAMIN_K
#وإلا
#define MAX_K	SAMIN_K
#endif
#وإلا
#ifdef XDOUBLE
#define MAX_K	XAMIN_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	ZAMIN_K
#وإلا
#define MAX_K	CAMIN_K
#endif
#endif

#endif

#وإلا

#define ABS

#ifndef USE_MIN

/* MAX */
#ifdef XDOUBLE
#define MAX_K	QMAX_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	DMAX_K
#وإلا
#define MAX_K	SMAX_K
#endif

#وإلا

/* MIN */
#ifdef XDOUBLE
#define MAX_K	QMIN_K
#وإلاإذا defined(DOUBLE)
#define MAX_K	DMIN_K
#وإلا
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

  𐤀𐤌 (n <= 0) أرجع 0;

#ifndef COMPLEX
  𐤀𐤌 (incx == 0) أرجع (ABS(*x));
#وإلا
  𐤀𐤌 (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  أرجع ret;
}

#وإلا

#ifdef COMPLEX
FLOAT CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  𐤀𐤌 (n <= 0) أرجع 0;

#ifndef COMPLEX
  𐤀𐤌 (incx == 0) أرجع (ABS(*x));
#وإلا
  𐤀𐤌 (incx == 0) أرجع (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  أرجع ret;
}

#endif
