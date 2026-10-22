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

#अगर_अघोषित USE_MIN

/* ABS & MAX */
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IQAMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IDAMAX_K
#अन्यथा
#घोषणा MAX_K	ISAMAX_K
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IXAMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IZAMAX_K
#अन्यथा
#घोषणा MAX_K	ICAMAX_K
#पूर्ण_अगर
#पूर्ण_अगर

#अन्यथा

/* ABS & MIN */
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IQAMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IDAMIN_K
#अन्यथा
#घोषणा MAX_K	ISAMIN_K
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IXAMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IZAMIN_K
#अन्यथा
#घोषणा MAX_K	ICAMIN_K
#पूर्ण_अगर
#पूर्ण_अगर

#पूर्ण_अगर

#अन्यथा

#अगर_अघोषित USE_MIN

/* MAX */
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IQMAX_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IDMAX_K
#अन्यथा
#घोषणा MAX_K	ISMAX_K
#पूर्ण_अगर

#अन्यथा

/* MIN */
#अगर_घोषित XDOUBLE
#घोषणा MAX_K	IQMIN_K
#या_अगर defined(DOUBLE)
#घोषणा MAX_K	IDMIN_K
#अन्यथा
#घोषणा MAX_K	ISMIN_K
#पूर्ण_अगर

#पूर्ण_अगर

#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS

blasint NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  blasint ret;

  PRINT_DEBUG_NAME;

  अगर (n <= 0) वापस 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (blasint)MAX_K(n, x, incx);

  अगर(ret > n) ret=n;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  वापस ret;
}

#अन्यथा
#अगर_घोषित COMPLEX
CBLAS_सूचकांक CNAME(blasint n, व्योम *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
CBLAS_सूचकांक CNAME(blasint n, FLOAT *x, blasint incx){
#पूर्ण_अगर

  CBLAS_सूचकांक ret;

  PRINT_DEBUG_CNAME;

  अगर (n <= 0) वापस 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  अगर (ret > n) ret=n;

  अगर (ret) ret --;

  FUNCTION_PROFILE_END(COMPSIZE, n, 0);

  IDEBUG_END;

  वापस ret;
}

#पूर्ण_अगर
