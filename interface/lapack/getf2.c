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

#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGETF2"
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGETF2"
#अन्यथा
#घोषणा ERROR_NAME "SGETF2"
#पूर्ण_अगर

OPENBLAS_EXPORT
व्योम NAME(blasint *M, blasint *N, FLOAT *a, blasint *ldA, blasint *ipiv, blasint *Info){

  blas_arg_t args;

  blasint info;
  FLOAT *बफर;
#अगर_घोषित PPC440
  बाह्य
#पूर्ण_अगर
  FLOAT *sa, *sb;

  PRINT_DEBUG_NAME;

  args.m    = *M;
  args.n    = *N;
  args.a    = (व्योम *)a;
  args.lda  = *ldA;
  args.c    = (व्योम *)ipiv;

  info  =    0;
  अगर (args.lda < MAX(1,args.m)) info = 4;
  अगर (args.n   < 0)             info = 2;
  अगर (args.m   < 0)             info = 1;
  अगर (info) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME) - 1);
    *Info = - info;
    वापस;
  }

  *Info = 0;
  अगर (args.m == 0 || args.n == 0) वापस;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

#अगर_अघोषित PPC440
  बफर = (FLOAT *)blas_memory_alloc(1);

  sa = (FLOAT *)((BLASLONG)बफर + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * आकार + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);
#पूर्ण_अगर

  info = GETF2(&args, शून्य, शून्य, sa, sb, 0);

  *Info = info;

#अगर_अघोषित PPC440
  blas_memory_free(बफर);
#पूर्ण_अगर

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.n,  2. / 3. * args.m * args.n * args.n);

  IDEBUG_END;

  वापस;
}
