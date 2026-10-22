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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#undef MAX_K

#ifdef USE_ABS

#જો defined(DOUBLE)
#define ABS fabs
#નહીંતર
#define ABS fabsf
#endif

#ifndef USE_MIN

/* ABS & MAX */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMAX_K
#અથવાજો defined(DOUBLE)
#define MAX_K	DAMAX_K
#નહીંતર
#define MAX_K	SAMAX_K
#endif
#નહીંતર
#ifdef XDOUBLE
#define MAX_K	XAMAX_K
#અથવાજો defined(DOUBLE)
#define MAX_K	ZAMAX_K
#નહીંતર
#define MAX_K	CAMAX_K
#endif
#endif

#નહીંતર

/* ABS & MIN */
#ifndef COMPLEX
#ifdef XDOUBLE
#define MAX_K	QAMIN_K
#અથવાજો defined(DOUBLE)
#define MAX_K	DAMIN_K
#નહીંતર
#define MAX_K	SAMIN_K
#endif
#નહીંતર
#ifdef XDOUBLE
#define MAX_K	XAMIN_K
#અથવાજો defined(DOUBLE)
#define MAX_K	ZAMIN_K
#નહીંતર
#define MAX_K	CAMIN_K
#endif
#endif

#endif

#નહીંતર

#define ABS

#ifndef USE_MIN

/* MAX */
#ifdef XDOUBLE
#define MAX_K	QMAX_K
#અથવાજો defined(DOUBLE)
#define MAX_K	DMAX_K
#નહીંતર
#define MAX_K	SMAX_K
#endif

#નહીંતર

/* MIN */
#ifdef XDOUBLE
#define MAX_K	QMIN_K
#અથવાજો defined(DOUBLE)
#define MAX_K	DMIN_K
#નહીંતર
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

  જો (n <= 0) પાછા 0;

#ifndef COMPLEX
  જો (incx == 0) પાછા (ABS(*x));
#નહીંતર
  જો (incx == 0) પાછા (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  પાછા ret;
}

#નહીંતર

#ifdef COMPLEX
FLOAT CNAME(blasint n, શૂન્ય *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#નહીંતર
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#endif
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  જો (n <= 0) પાછા 0;

#ifndef COMPLEX
  જો (incx == 0) પાછા (ABS(*x));
#નહીંતર
  જો (incx == 0) પાછા (ABS(*x) + ABS(*(x+1)));
#endif

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  પાછા ret;
}

#endif
