#ifdef CBLAS

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <stdlib.h>
#অন্তর্ভাব <string.h>
#অন্তর্ভাব <stdarg.h>
#অন্তর্ভাব "common.h"

OPENBLAS_EXPORT

শূন্য CNAME(blasint p, বর্ণ *rout, বর্ণ *form, ...)
{
   va_list args;

   va_start(args, form);

   যদি (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

