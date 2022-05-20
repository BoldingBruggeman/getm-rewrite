Here, the budget equation for layer-averaged momentum in northern direction,
$q_k$,
is calculated. The physical equation is given as equation (\ref{vEq}),
the layer-integrated equation as (\ref{vEqvi}), and after curvilinear
transformation as (\ref{vEqviCurvi}).
In this routine, first the Coriolis rotation term, $fp_k$ is calculated,
either as direct transport averaging, or following \cite{ESPELIDea00}
by using velocity averages (in case the compiler option {\tt NEW\_CORI}
is set).

As a next step, explicit forcing terms (advection, diffusion,
internal pressure gradient, surface stresses) are added up (into the variable
{\tt ex(k)}), the eddy viscosity is horizontally interpolated to the V-point,
and the barotropic pressure gradient is calculated (the latter
includes the pressure gradient correction for drying points, see
section \ref{Section_dry}).
Afterwards, the matrix is set up for each water column, and it is solved
by means of a tri-diagonal matrix solver.

In case that the compiler option {\tt STRUCTURE\_FRICTION} is switched on,
the frictional effect of structures in the water column is calculated
by adding the quadratic frictional term $C v \sqrt{u^2+v^2}$ (with a minus
sign on
the right hand side) numerically implicitly to the $v$-equation,
with the friction coefficient $C$. The explicit part of this term,
$C\sqrt{u^2+v^2}$,
is calculated in the routine {\tt structure\_friction\_3d.F90}.

Finally, the new velocity profile is shifted such that its vertical
integral is identical to the time integral of the vertically integrated
transport.
If the compiler option {\tt MUDFLAT} is defined, this fitting of profiles
is made with
respect to the new surface elevation, otherwise to the
old surface elevation.
