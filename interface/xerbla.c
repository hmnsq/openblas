#ifdef CBLAS

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ <stdlib.h>
#অন্তর্ভাৱ <string.h>
#অন্তর্ভাৱ <stdarg.h>
#অন্তর্ভাৱ "common.h"

OPENBLAS_EXPORT

শূন্য CNAME(blasint p, ৱর্ণ *rout, ৱর্ণ *form, ...)
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

