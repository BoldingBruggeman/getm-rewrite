module advection_superbee

#define SUPERBEE
#define _LIMITER1_ limiter = MAX( MIN( 2*ratio , 1.0_real64 ) , MIN( ratio , 2.0_real64 ) )
#define _LIMITER2_
#define _LIMITER3_
#define _LIMITER4_
#define _TYPE_NAME_ type_advection_superbee

#include "advection.F90.template"

end module
