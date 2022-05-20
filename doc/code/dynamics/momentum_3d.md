### Three-dimensional momentum equations

For geophysical coastal sea
and ocean dynamics, usually the three-dimensional hydrostatic equations
of motion with the Boussinesq approximation and the eddy
viscosity assumption are used (\cite{BRYAN69}, \cite{COX84},
\cite{BLUMBERGea87}, \cite{HAIDVOGELea99}, \cite{KANTHAea00b}).
In the flux form, the dynamic equations of motion for
the horizontal velocity components can be written in Cartesian
coordinates as:


\begin{equation}\label{uEq}
\begin{array}{l}
\displaystyle
\partial_t u
+\partial_z(uw)
-\partial_z\left((\nu_t+\nu) \partial_z u\right)
\\ \\ \displaystyle
\qquad+\alpha\bigg(\partial_x(u^2)+\partial_y(uv)
-\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)
\\ \\ \displaystyle
\qquad
-fv
-\int_z^{\zeta}\partial_x b\,dz' \bigg)
=
- g\partial_x \zeta,
\end{array}
\end{equation}

\begin{equation}\label{vEq}
\begin{array}{l}
\displaystyle
\partial_t v +\partial_z(vw)
-\partial_z\left((\nu_t+\nu) \partial_z v\right)
\\ \\
\displaystyle
\qquad+\alpha\bigg(\partial_x(vu)+\partial_y(v^2)
-\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)
\\ \\
\displaystyle
\qquad
+fu
-\int_z^{\zeta}\partial_y b\,dz' \bigg)=
- g\partial_y \zeta.
\end{array}
\end{equation}

The vertical velocity is calculated by means of the
incompressibility condition:

\begin{equation}\label{Konti}
\partial_x u +\partial_y v +\partial_z w = 0.
\end{equation}


Here, $u$, $v$ and $w$ are the ensemble averaged
velocity components with respect
to the $x$, $y$ and $z$ direction, respectively.
The vertical coordinate $z$ ranges from the bottom $-H(x,y)$ \label{Hxy}
to the surface $\zeta(t,x,y)$ with $t$ denoting time.
$\nu_t$ is the vertical eddy viscosity, $\nu$ the kinematic viscosity,
$f$ the Coriolis
parameter, and $g$ is the gravitational acceleration.
The horizontal mixing is parameterised by terms containing the
horizontal eddy viscosity $A_h^M$, see \cite{BLUMBERGea87}.
The buoyancy $b$ is defined as

\begin{equation}\label{bdef}
b=-g\frac{\rho-\rho_0}{\rho_0}
\end{equation}

with the density $\rho$ and a reference density $\rho_0$.
The last term on the left hand sides of equations (\ref{uEq}) and (\ref{vEq})
are the internal (due to density gradients)
and the terms on the right hand sides are the external
(due to surface slopes) pressure gradients. In the latter, the deviation of
surface density from reference density is neglected (see \cite{BURCHARDea97}).
The derivation of equations (\ref{uEq}) - (\ref{Konti}) has been shown in
numerous publications, see e.g.\ \cite{PEDLOSKY87}, \cite{HAIDVOGELea99},
\cite{BURCHARD02}.

In hydrostatic 3D models, the vertical velocity is calculated by means of
equation (\ref{Konti}) velocity equation.
Due to this, mass conservation and free surface elevation
can easily be obtained.

Drying and flooding of mud-flats is already incorporated in
the physical equations by multiplying some terms with the
non-dimensional number $\alpha$ which equals unity in regions where a
critical water depth $D_{crit}$ is exceeded and approaches zero
when the water depth $D$ tends to a minimum value $D_{min}$:

\begin{equation}\label{alpha}
\alpha=\min\left\{1,\frac{D-D_{min}}{D_{crit}-D_{min}}\right\}.
\end{equation}

Thus, $\alpha=1$ for $D\geq D_{crit}$,  such that the usual momentum
equation results except for very shallow water, where simplified physics
are considered with a balance between tendency, friction and external pressure
gradient. In a typical wadden sea application, $D_{crit}$ is of the order
of 0.1 m and $D_{\min}$ of the order of 0.02 m (see \cite{BURCHARD98},
\cite{BURCHARDea03a}).

### Kinematic boundary conditions and surface elevation equation

At the surface and at the bottom, kinematic boundary conditions
result from the requirement that
the particles at the boundaries are moving along these boundaries:

\begin{equation}\label{KinBoundsSurf}
w= \partial_t \zeta +  u \partial_x\zeta + v \partial_y \zeta
\qquad \mbox{for } z=\zeta,
\end{equation}

\begin{equation}\label{KinBoundsBott}
w= - u \partial_xH - v \partial_y H
\qquad \mbox{for } z=-H.
\end{equation}

### Dynamic boundary conditions

At the bottom boundaries, no-slip conditions are prescribed for the horizontal
velocity components:

\begin{equation}\label{DynBoundbx}
u=0, \quad v=0.
\end{equation}

With (\ref{KinBoundsBott}), also $w=0$ holds at the bottom.
It should be noted already here, that the bottom boundary condition
(\ref{DynBoundbx}) is generally not directly used in numerical ocean models,
since the near-bottom values of the horizontal velocity components
are not located at the bed, but half a grid box above it.
Instead, a logarithmic velocity profile is assumed in the bottom layer,
leading to a quadratic friction law, see section \ref{sec-bottom-friction-3d}.

At the surface, the dynamic boundary conditions read:

\begin{equation}\label{DynBoundSx}
\begin{array}{l}
(\nu_t+\nu) \partial_z u=\alpha \tau^x_{s},
\\ \\
(\nu_t+\nu) \partial_z v=\alpha \tau^y_{s},
\end{array}
\end{equation}

The surface stresses (normalised by the reference density) $\tau_s^x$
and $\tau_s^y$ are calculated as functions of wind speed, wind direction,
surface roughness etc.
Also here, the drying parameter $\alpha$ is included in order to
provide an easy handling of drying and flooding.

### Lateral boundary conditions

Let $G$\label{G} denote the lateral boundary of the model domain with the
closed land boundary $G^c$\label{Gc} and the open boundary $G^o$\label{Go} such that
$G^c \cup G^o=G$ and $G^c \cap G^o=\emptyset$.
Let further $\vec u=(u,v)$\label{vecu} denote the horizontal velocity vector
and $\vec u_n=(-v,u)$\label{vecun} its normal vector.
At closed boundaries, the flow must be parallel to the boundary:

\begin{equation}
\vec u_n \cdot \vec\nabla G^c = 0 
\end{equation}

with $\vec\nabla=(\partial_x,\partial_y)$ being the gradient operator.

For an eastern or a western closed boundary with
$\vec\nabla G^c=(0,1)$ this has the consequence that
$u=0$ and, equivalently, for a southern
or a northern closed boundary with $\vec\nabla G^c=(1,0)$
this has the consequence that $v=0$.

At open boundaries, the velocity gradients across the boundary vanish:

\begin{equation}
\vec\nabla_n u \cdot \vec\nabla G^o = 0, \quad 
\vec\nabla_n v \cdot \vec\nabla G^o = 0, \quad 
\end{equation}

with $\vec\nabla_n=(-\partial_y,\partial_x)$ being the
operator normal to the gradient operator.

For an eastern or a western open boundary with
this has the consequence that
$\partial_x u=\partial_x v =0$ and, equivalently, for a southern
or a northern open boundary
this has the consequence that $\partial_y u=\partial_y v =0$.

At so-called forced
open boundaries, the sea surface elevation $\zeta$ is prescribed.
At passive open boundaries, it is assumed that the curvature of the
surface elevation normal to the boundary is zero, with the consequence that
the spatial derivatives of the surface slopes normal to the boundaries vanish.

### Layer-integrated equations

There are two different ways to derive the layer-integrated
equations. \cite{BURCHARDea97} transform first the
equations into general vertical coordinate form
(see \cite{DELEERSNIJDERea92})
and afterwards
integrate the transformed equations over constant intervals in
the transformed space. \cite{LANDERea94} integrate the
equations in the Cartesian space over surfaces $z_k$ by
considering the Leibniz rule

\begin{equation}\label{Leibniz}
\int_{z_{k-1}}^{z_k}\partial_x f\,dz
=
\partial_x\int_{z_{k-1}}^{z_k} f\,dz
-f(z_k)\partial_xz_k
+f(z_{k-1})\partial_xz_{k-1}
\end{equation}

for any function $f$.
For the vertical staggering of the layer notation see figure
\ref{figvertgrid}.

More details about the layer integration are given in
\cite{BURCHARDea97}.

With the further definitions of layer integrated transport,

\begin{equation}\label{pqdef} 
p_k:=\int_{z_{k-1}}^{z_k}u\,dz,\qquad
q_k:=\int_{z_{k-1}}^{z_k}v\,dz,
\end{equation}

layer mean velocities,

\begin{equation}\label{ukvkdef} 
u_k:=\frac{p_k}{h_k},\qquad v_k:=\frac{q_k}{h_k},
\end{equation}

and layer averaged tracer concentrations and buoyancy,

\begin{equation}\label{ckbkdef} 
c^i_k:=\frac{1}{h_k}\int_{z_{k-1}}^{z_k}c^i\,dz,\qquad
b_k:=\frac{1}{h_k}\int_{z_{k-1}}^{z_k}b\,dz,
\end{equation}

and the grid related vertical velocity,

\begin{equation}\label{barwdef} 
\bar w_k:=(w-\partial_tz-u\partial_xz-v\partial_yz)_{z=z_k},
\end{equation}

the continuity equation (\ref{Konti}) has the layer-integrated form:

\begin{equation}\label{ContiLayerInt}
\partial_t h_k + \partial_x p_k + \partial_y q_k + \bar w_k - \bar w_{k-1}=0.
\end{equation}

It should be noted that the grid related velocity is located on the
layer interfaces.
After this, the layer-integrated momentum equations read as:

\begin{equation}\label{uEqvi}
\begin{array}{l}
\partial_t p_k 
+\bar w_k \tilde u_k -\bar w_{k-1} \tilde u_{k-1} 
-\tau^x_k + \tau^x_{k-1} 
\\ \\ \quad
+\alpha\Bigg\{\partial_x(u_kp_k)+\partial_y(v_kp_k)
\\ \\ \displaystyle \quad 
-\partial_x\left(2A_k^Mh_k\partial_xu_k\right)-\partial_y\left(A_k^Mh_k
(\partial_yu_k+\partial_xv_k)\right)
- fq_k 
\\ \\ \quad
\displaystyle
-h_k\left(
\frac12h_N(\partial^*_xb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_xb)_j
\right)\Bigg\}
=
-gh_k\partial_x\zeta,
\end{array}
\end{equation}

\begin{equation}\label{vEqvi}
\begin{array}{l}
\partial_t q_k 
+\bar w_k  \tilde v_k -\bar w_{k-1}  \tilde v_{k-1} 
-\tau^y_k + \tau^y_{k-1} 
\\ \\ \quad
+\alpha\Bigg\{\partial_x(u_kq_k)+\partial_y(v_kq_k)
\\ \\ \displaystyle \quad 
-\partial_y\left(2A_k^Mh_k\partial_yv_k\right)-\partial_x\left(A_k^Mh_k
(\partial_yu_k+\partial_xv_k)\right)
+ fp_k 
\\ \\ \quad
\displaystyle
-h_k\left(
\frac12h_N(\partial^*_yb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_yb)_j
\right)\Bigg\}
=
-gh_k\partial_y\zeta
\end{array}
\end{equation}

with suitably chosen advective horizontal
velocities $\tilde u_k$ and $\tilde v_k$
(see section \ref{sec-uv-advect-3d}) on page \pageref{sec-uv-advect-3d},
the shear stresses

\begin{equation}\label{tauxkdef}
\tau^x_k = \left(\nu_t \partial_z u \right)_k,
\end{equation}

and

\begin{equation}\label{tauykdef}
\tau^y_k = \left(\nu_t \partial_z v \right)_k,
\end{equation}

and the horizontal buoyancy gradients

\begin{equation}
(\partial^*_xb)_k=\frac12(\partial_xb_{k+1}+\partial_x b_k)
-\partial_xz_k\frac{b_{k+1}-b_k}{\frac12(h_{k+1}+h_k)}
\end{equation}

and

\begin{equation}
(\partial^*_yb)_k=\frac12(\partial_yb_{k+1}+\partial_y b_k)
-\partial_yz_k\frac{b_{k+1}-b_k}{\frac12(h_{k+1}+h_k)}.
\end{equation}

The layer integration of the pressure gradient force
is discussed in detail by \cite{BURCHARDea97}.

A conservative formulation can be derived
for the recalculation of the physical vertical velocity
$w$ which is convenient
in the discrete space if $w$ is evaluated at the layer centres
(see \cite{DELEERSNIJDERea92}):

\begin{equation}\label{conservative_w}
w_k=\frac{1}{h_k}\left(
\partial_t(h_kz_{k-1/2})+\partial_x(p_kz_{k-1/2})+\partial_y(q_kz_{k-1/2})
+\bar w_kz_k-\bar w_{k-1}z_{k-1}\right).
\end{equation}

It should be mentioned that $w$ only needs to be evaluated for
post-processing reasons.

For the layer-integrated tracer concentrations, we obtain
the following expression:

\begin{equation}\label{C_Layer_Int}
\begin{array}{l}
\partial_t (h_k c^i_k) + \partial_x (p_kc^i_k)+\partial_y (q_k c^i_k)
+(\bar w_k+w^s_k) \tilde c^i_k 
- (\bar w_{k-1}+w^s_{k-1}) \tilde c^i_{k-1}\\
\\ \qquad
-(\nu'_t\partial_z c^i)_{k}
+(\nu'_t\partial_z c^i)_{k-1}
-\partial_x\left(A_k^Th_k\partial_xc^i_k\right)
-\partial_y\left(A_k^Th_k\partial_yc^i_k\right)
=Q^i_k.
\end{array}
\end{equation}

It should be noted that the "horizontal" diffusion does no longer occur
along geopotential surfaces but along horizontal coordinate lines.
The properly transformed formulation would include some cross-diagonal
terms which may lead to numerical instabilities due to violation
of monotonicity. For an in-depth discussion of this
problem, see \cite{BECKERSea98} and \cite{BECKERSea00}.

