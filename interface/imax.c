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

#यदि_अपरिभाषितम् USE_MIN

/* ABS & MAX */
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IQAMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IDAMAX_K
#अन्यथा
#परिभाषय MAX_K	ISAMAX_K
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IXAMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IZAMAX_K
#अन्यथा
#परिभाषय MAX_K	ICAMAX_K
#यदि_अन्तः
#यदि_अन्तः

#अन्यथा

/* ABS & MIN */
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IQAMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IDAMIN_K
#अन्यथा
#परिभाषय MAX_K	ISAMIN_K
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IXAMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IZAMIN_K
#अन्यथा
#परिभाषय MAX_K	ICAMIN_K
#यदि_अन्तः
#यदि_अन्तः

#यदि_अन्तः

#अन्यथा

#यदि_अपरिभाषितम् USE_MIN

/* MAX */
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IQMAX_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IDMAX_K
#अन्यथा
#परिभाषय MAX_K	ISMAX_K
#यदि_अन्तः

#अन्यथा

/* MIN */
#यदि_परिभाषितम् XDOUBLE
#परिभाषय MAX_K	IQMIN_K
#नोचेद्यदि defined(DOUBLE)
#परिभाषय MAX_K	IDMIN_K
#अन्यथा
#परिभाषय MAX_K	ISMIN_K
#यदि_अन्तः

#यदि_अन्तः

#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

blasint NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  blasint ret;

  PRINT_DEBUG_NAME;

  यदि (n <= 0) प्रत्यागम 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = (blasint)MAX_K(n, x, incx);

  यदि(ret > n) ret=n;

  FUNCTION_PROFILE_अन्तः(COMPSIZE, n, 0);

  IDEBUG_अन्तः;

  प्रत्यागम ret;
}

#अन्यथा
#यदि_परिभाषितम् COMPLEX
CBLAS_सूचकांक CNAME(blasint n, शून्य *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#अन्यथा
CBLAS_सूचकांक CNAME(blasint n, FLOAT *x, blasint incx){
#यदि_अन्तः

  CBLAS_सूचकांक ret;

  PRINT_DEBUG_CNAME;

  यदि (n <= 0) प्रत्यागम 0;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  ret = MAX_K(n, x, incx);

  यदि (ret > n) ret=n;

  यदि (ret) ret --;

  FUNCTION_PROFILE_अन्तः(COMPSIZE, n, 0);

  IDEBUG_अन्तः;

  प्रत्यागम ret;
}

#यदि_अन्तः
