#define MUSCL

#define _MODULE_NAME_ advection_muscl
#define _U_SUB_NAME_ u_advection_muscl
#define _V_SUB_NAME_ v_advection_muscl
#define _W_SUB_NAME_ w_advection_muscl

#include "advection.F90.template"
#undef MUSCL
