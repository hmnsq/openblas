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
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

#अपरिभाषय MAX_K

#यदि_परिभाषितम् USE_ABS

#यदि defined(DOUBLE)
#परिभाषय ABS fabs
#अन्यथा
#परिभाषय ABS fabsf
#यदि_अन्तः

#यदि_अपरिभाषितम् USE_MIN

/* ABS & MAX */
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	QAMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	DAMAX_K
#अन्यथा
#परिभाषय MAX_K	SAMAX_K
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	XAMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	ZAMAX_K
#अन्यथा
#परिभाषय MAX_K	CAMAX_K
#यदि_अन्तः
#यदि_अन्तः

#अन्यथा

/* ABS & MIN */
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	QAMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	DAMIN_K
#अन्यथा
#परिभाषय MAX_K	SAMIN_K
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	XAMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	ZAMIN_K
#अन्यथा
#परिभाषय MAX_K	CAMIN_K
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः

#अन्यथा

#परिभाषय ABS

#यदि_अपरिभाषितम् USE_MIN

/* MAX */
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	QMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	DMAX_K
#अन्यथा
#परिभाषय MAX_K	SMAX_K
#यदि_अन्तः

#अन्यथा

/* MIN */
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	QMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	DMIN_K
#अन्यथा
#परिभाषय MAX_K	SMIN_K
#यदि_अन्तः

#यदि_अन्तः

#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  यदि (n <= 0) प्रत्यागम 0;

#यदि_अपरिभाषितम् COMPLEX
  यदि (incx == 0) प्रत्यागम (ABS(*x));
#अन्यथा
  यदि (incx == 0) प्रत्यागम (ABS(*x) + ABS(*(x+1)));
#यदि_अन्तः

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_अन्तः(COMPSIZE, n, 0);

  IDEBUG_अन्तः;

  प्रत्यागम ret;
}

#अन्यथा

#यदि_परिभाषितम् COMPLEX
FLOAT CNAME(blasint n, शून्य *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#यदि_अन्तः
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  यदि (n <= 0) प्रत्यागम 0;

#यदि_अपरिभाषितम् COMPLEX
  यदि (incx == 0) प्रत्यागम (ABS(*x));
#अन्यथा
  यदि (incx == 0) प्रत्यागम (ABS(*x) + ABS(*(x+1)));
#यदि_अन्तः

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_अन्तः(COMPSIZE, n, 0);

  IDEBUG_अन्तः;

  प्रत्यागम ret;
}

#यदि_अन्तः
