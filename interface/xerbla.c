#ifdef CBLAS

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ <string.h>
#അന്തര്ഭാവ <stdarg.h>
#അന്തര്ഭാവ "common.h"

OPENBLAS_EXPORT

ശൂന്യം CNAME(blasint p, അക്ഷരം *rout, അക്ഷരം *form, ...)
{
   va_list args;

   va_start(args, form);

   എങ്കിൽ (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

