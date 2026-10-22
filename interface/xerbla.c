#ifdef CBLAS

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ <string.h>
#ଅନ୍ତର୍ଭାଵ <stdarg.h>
#ଅନ୍ତର୍ଭାଵ "common.h"

OPENBLAS_EXPORT

ଶୂନ୍ୟ CNAME(blasint p, ଅକ୍ଷର *rout, ଅକ୍ଷର *form, ...)
{
   va_list args;

   va_start(args, form);

   ଯଦି (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

