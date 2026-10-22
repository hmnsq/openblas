#ifdef CBLAS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "common.h"

OPENBLAS_EXPORT

void CNAME(blasint p, char *rout, char *form, ...)
{
   va_list मगा_कातविरान;

   va_start(मगा_कातविरान, form);

   if (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, मगा_कातविरान);
   va_end(मगा_कातविरान);
   exit(-1);
}
#endif

