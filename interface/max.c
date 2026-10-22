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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਪਰਿਭਾਸ਼ਾ_ਰੱਦ MAX_K

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ USE_ABS

#ਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ABS fabs
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ABS fabsf
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ USE_MIN

/* ABS & MAX */
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	QAMAX_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	DAMAX_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	SAMAX_K
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	XAMAX_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	ZAMAX_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	CAMAX_K
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਨਹੀਂਤਾਂ

/* ABS & MIN */
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	QAMIN_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	DAMIN_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	SAMIN_K
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	XAMIN_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	ZAMIN_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	CAMIN_K
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ

#ਨਹੀਂਤਾਂ

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ABS

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ USE_MIN

/* MAX */
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	QMAX_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	DMAX_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	SMAX_K
#ਜੇ_ਅੰਤ

#ਨਹੀਂਤਾਂ

/* MIN */
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	QMIN_K
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	DMIN_K
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ MAX_K	SMIN_K
#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  ਜੇ (n <= 0) ਵਾਪਸ 0;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (incx == 0) ਵਾਪਸ (ABS(*x));
#ਨਹੀਂਤਾਂ
  ਜੇ (incx == 0) ਵਾਪਸ (ABS(*x) + ABS(*(x+1)));
#ਜੇ_ਅੰਤ

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE, n, 0);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;
}

#ਨਹੀਂਤਾਂ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX
FLOAT CNAME(blasint n, ਖਾਲੀ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#ਨਹੀਂਤਾਂ
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#ਜੇ_ਅੰਤ
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  ਜੇ (n <= 0) ਵਾਪਸ 0;

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (incx == 0) ਵਾਪਸ (ABS(*x));
#ਨਹੀਂਤਾਂ
  ਜੇ (incx == 0) ਵਾਪਸ (ABS(*x) + ABS(*(x+1)));
#ਜੇ_ਅੰਤ

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE, n, 0);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ ret;
}

#ਜੇ_ਅੰਤ
