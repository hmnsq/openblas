#ifdef CBLAS

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <stdlib.h>
#அந்தர்பாவ <string.h>
#அந்தர்பாவ <stdarg.h>
#அந்தர்பாவ "common.h"

OPENBLAS_EXPORT

வெறுமை CNAME(blasint p, எழுத்து *rout, எழுத்து *form, ...)
{
   va_list args;

   va_start(args, form);

   எனில் (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

