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

#समावेश <stdio.h>
#समावेश "common.h"
#अगर_घोषित FUNCTION_PROFILE
#समावेश "functable.h"
#पूर्ण_अगर

#अघोषित MAX_K

#अगर_घोषित USE_ABS

#अगर defined(DOUBLE)
#घोषणा ABS भ_असल
#अन्यथा
#घोषणा ABS fabsf
#पूर्ण_अगर

#अगर_अघोषित USE_MIN

/* ABS & MAX */
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	QAMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	DAMAX_K
#अन्यथा
#घोषणा MAX_K	SAMAX_K
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	XAMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	ZAMAX_K
#अन्यथा
#घोषणा MAX_K	CAMAX_K
#पूर्ण_अगर
#पूर्ण_अगर

#अन्यथा

/* ABS & MIN */
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	QAMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	DAMIN_K
#अन्यथा
#घोषणा MAX_K	SAMIN_K
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	XAMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	ZAMIN_K
#अन्यथा
#घोषणा MAX_K	CAMIN_K
#पूर्ण_अगर
#पूर्ण_अगर

#पूर्ण_अगर

#अन्यथा

#घोषणा ABS

#अगर_अघोषित USE_MIN

/* MAX */
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	QMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	DMAX_K
#अन्यथा
#घोषणा MAX_K	SMAX_K
#पूर्ण_अगर

#अन्यथा

/* MIN */
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	QMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	DMIN_K
#अन्यथा
#घोषणा MAX_K	SMIN_K
#पूर्ण_अगर

#पूर्ण_अगर

#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  FLOATRET ret;

  PRINT_DEBUG_NAME;

  अगर (n <= 0) वापस 0;

#अगर_अघोषित COMPLEX
  अगर (incx == 0) वापस (ABS(*x));
#अन्यथा
  अगर (incx == 0) वापस (ABS(*x) + ABS(*(x+1)));
#पूर्ण_अगर

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (FLOATRET)MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  वापस ret;
}

#अन्यथा

#अगर_घोषित COMPLEX
FLOAT CNAME(blasint n, व्योम *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
FLOAT CNAME(blasint n, FLOAT *x, blasint incx){
#पूर्ण_अगर
  
  FLOAT ret;

  PRINT_DEBUG_CNAME;

  अगर (n <= 0) वापस 0;

#अगर_अघोषित COMPLEX
  अगर (incx == 0) वापस (ABS(*x));
#अन्यथा
  अगर (incx == 0) वापस (ABS(*x) + ABS(*(x+1)));
#पूर्ण_अगर

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  वापस ret;
}

#पूर्ण_अगर
