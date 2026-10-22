#ifdef CBLAS

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <stdlib.h>
#અન્તર્ભાવ <string.h>
#અન્તર્ભાવ <stdarg.h>
#અન્તર્ભાવ "common.h"

OPENBLAS_EXPORT

શૂન્ય CNAME(blasint p, અક્ષર *rout, અક્ષર *form, ...)
{
   va_list args;

   va_start(args, form);

   જો (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

