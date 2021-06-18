module advection_upstream

#define  UPSTREAM
#define _LIMITER1_ limiter = 0._real64
#define _LIMITER2_
#define _LIMITER3_
#define _LIMITER4_
#define _TYPE_NAME_ type_advection_upstream

#include "advection.F90.template"

end module
