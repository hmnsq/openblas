#నిర్వచితమైతే CBLAS

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <stdlib.h>
#అన్తర్భావ <string.h>
#అన్తర్భావ <stdarg.h>
#అన్తర్భావ "common.h"

OPENBLAS_EXPORT

శూన్య CNAME(blasint p, అక్షరం *rout, అక్షరం *form, ...)
{
   va_list వాదనలు;

   va_start(వాదనలు, form);

   ఒకవేళ (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (దీర్ఘ దీర్ఘ)p, rout);
   vfprintf(stderr, form, వాదనలు);
   va_end(వాదనలు);
   నిష్క్రమించు(-1);
}
#అయితే_ముగింపు

