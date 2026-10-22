#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ CBLAS

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdlib.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <string.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdarg.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"

OPENBLAS_EXPORT

ਖਾਲੀ CNAME(blasint p, ਅੱਖਰ *rout, ਅੱਖਰ *form, ...)
{
   va_list ਦਲੀਲਾਂ;

   va_start(ਦਲੀਲਾਂ, form);

   ਜੇ (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (ਲੰਮਾ ਲੰਮਾ)p, rout);
   vfprintf(stderr, form, ਦਲੀਲਾਂ);
   va_end(ਦਲੀਲਾਂ);
   ਬਾਹਰ_ਨਿਕਲੋ(-1);
}
#ਜੇ_ਅੰਤ

