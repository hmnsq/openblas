#ifdef CBLAS

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل <string.h>
#شامل <stdarg.h>
#شامل "common.h"

OPENBLAS_EXPORT

تهی CNAME(blasint p, نویسه *rout, نویسه *form, ...)
{
   va_list args;

   va_start(args, form);

   اگر (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

