#define SUPERBEE
#define _LIMITER1_ limiter = MAX( MIN( 2*ratio , 1.0_real64 ) , MIN( ratio , 2.0_real64 ) )
#define _LIMITER2_
#define _LIMITER3_
#define _LIMITER4_

#define _MODULE_NAME_ advection_superbee
#define _U_SUB_NAME_ u_advection_superbee
#define _V_SUB_NAME_ v_advection_superbee
#define _W_SUB_NAME_ W_advection_superbee

#include "advection.F90.template"
