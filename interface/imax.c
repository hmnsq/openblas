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
#که_تعریف_شوی FUNCTION_PROFILE
#ضمن "functable.h"
#د_که_پای

#تعریف_لغوه MAX_K

#که_تعریف_شوی USE_ABS

#که_نه_تعریف_شوی USE_MIN

/* ABS & MAX */
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IQAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IDAMAX_K
#وإلا
#تعریف_کړه MAX_K	ISAMAX_K
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IXAMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IZAMAX_K
#وإلا
#تعریف_کړه MAX_K	ICAMAX_K
#د_که_پای
#د_که_پای

#وإلا

/* ABS & MIN */
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IQAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IDAMIN_K
#وإلا
#تعریف_کړه MAX_K	ISAMIN_K
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IXAMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IZAMIN_K
#وإلا
#تعریف_کړه MAX_K	ICAMIN_K
#د_که_پای
#د_که_پای

#د_که_پای

#وإلا

#که_نه_تعریف_شوی USE_MIN

/* MAX */
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IQMAX_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IDMAX_K
#وإلا
#تعریف_کړه MAX_K	ISMAX_K
#د_که_پای

#وإلا

/* MIN */
#که_تعریف_شوی XDOUBLE
#تعریف_کړه MAX_K	IQMIN_K
#وإلاإذا defined(DOUBLE)
#تعریف_کړه MAX_K	IDMIN_K
#وإلا
#تعریف_کړه MAX_K	ISMIN_K
#د_که_پای

#د_که_پای

#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS

blasint NAME(blasint *N, FLOAT *x, blasint *INCX){

  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  blasint ret;

  PRINT_DEBUG_NAME;

  إذا (n <= 0) أرجع 0;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  ret = (blasint)MAX_K(n, x, incx);

  إذا(ret > n) ret=n;

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE, n, 0);

  IDEBUG_پای_ټکی;

  أرجع ret;
}

#وإلا
#که_تعریف_شوی COMPLEX
CBLAS_شاخص CNAME(blasint n, فراغ *vx, blasint incx){
  FLOAT *x = (FLOAT*) vx;
#وإلا
CBLAS_شاخص CNAME(blasint n, FLOAT *x, blasint incx){
#د_که_پای

  CBLAS_شاخص ret;

  PRINT_DEBUG_CNAME;

  إذا (n <= 0) أرجع 0;

  IDEBUG_پیل;

  FUNCTION_PROFILE_پیل();

  ret = MAX_K(n, x, incx);

  إذا (ret > n) ret=n;

  إذا (ret) ret --;

  FUNCTION_PROFILE_پای_ټکی(COMPSIZE, n, 0);

  IDEBUG_پای_ټکی;

  أرجع ret;
}

#د_که_پای
