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

#ifdef XDOUBLE
#define ERROR_NAME "QTRTI2"
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DTRTI2"
#अन्यथा
#define ERROR_NAME "STRTI2"
#endif

स्थिर blasint (*trti2[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifdef XDOUBLE
  qtrti2_UU, qtrti2_UN, qtrti2_LU, qtrti2_LN,
#किंवाजर defined(DOUBLE)
  dtrti2_UU, dtrti2_UN, dtrti2_LU, dtrti2_LN,
#अन्यथा
  strti2_UU, strti2_UN, strti2_LU, strti2_LN,
#endif
  };

OPENBLAS_EXPORT
शून्य NAME(अक्षर *UPLO, अक्षर *DIAG, blasint *N, FLOAT *a, blasint *ldA, blasint *Info){

  blas_arg_t args;

  blasint uplo_arg = *UPLO;
  blasint diag_arg = *DIAG;
  blasint uplo, diag;
  blasint info;
  FLOAT *buffer;
#ifdef PPC440
  बाह्य
#endif
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.n    = *N;
  args.a    = (शून्य *)a;
  args.lda  = *ldA;

  TOUPPER(uplo_arg);
  TOUPPER(diag_arg);

  uplo = -1;
  जर (uplo_arg == 'U') uplo = 0;
  जर (uplo_arg == 'L') uplo = 1;
  diag = -1;
  जर (diag_arg == 'U') diag = 0;
  जर (diag_arg == 'N') diag = 1;

  info  = 0;
  जर (args.lda  < MAX(1,args.n)) info = 5;
  जर (args.n    < 0)             info = 3;
  जर (diag < 0)                  info = 2;
  जर (uplo < 0)                  info = 1;
  जर (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME) - 1);
    *Info = - info;
    परत;
  }

  *Info = 0;

  जर (args.n <= 0) परत;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#ifndef PPC440
  buffer = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#endif

  info = (trti2[(uplo << 1) | diag])(&args, रिक्त, रिक्त, sa, sb, 0);

  *Info = info;

#ifndef PPC440
  blas_memory_free(buffer);
#endif

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, .5 * args.n * args.n,
		          args.n * (1./3. + args.n * ( 1./2. + args.n * 1./6.))
		       +  args.n * (1./3. + args.n * (-1./2. + args.n * 1./6.)));

  IDEBUG_END;

  परत;
}
