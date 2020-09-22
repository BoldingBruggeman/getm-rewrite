#define  UPSTREAM
#define _LIMITER1_ limiter = 0._real64
#define _LIMITER2_
#define _LIMITER3_
#define _LIMITER4_

#define _MODULE_NAME_ advection_upstream
#define _U_SUB_NAME_ u_advection_upstream
#define _V_SUB_NAME_ v_advection_upstream
#define _W_SUB_NAME_ w_advection_upstream

#include "advection.F90.template"
#undef UPSTREAM
