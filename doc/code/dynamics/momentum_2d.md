## Vertically integrated mode

In order to provide the splitting of the model into
an internal and an external mode, the continuity equation and
the momentum equations are vertically integrated.
The vertical integral of the continuity equation together
with the kinematic boundary conditions (\ref{KinBoundsSurf}) and
(\ref{KinBoundsBott}) gives the sea surface elevation equation:

\begin{equation}\label{Elevation}
\partial_t \zeta = - \partial_x U- \partial_y V.
\end{equation}

with

\begin{equation}
U=\int_{-H}^{\zeta} u\,dz,\qquad V=\int_{-H}^{\zeta} v\,dz.
\end{equation}

Integrating the momentum equations (\ref{uEq}) and (\ref{vEq}) vertically
results in:

\begin{equation}\label{Uint}
\begin{array}{l}
\displaystyle
\partial_tU+
\tau_b^x+\alpha\bigg(\int_{-H}^{\zeta}
\left(\partial_x u^2+\partial_y(uv)\right)\,dz\\
\\
\displaystyle
\quad-\tau^x_s
-\int_{-H}^{\zeta}\big(
\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle
\quad -fV
-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_x b\,dz'\,dz \bigg)
=
-gD\partial_x\zeta
\end{array}
\end{equation}

and


\begin{equation}\label{Vint}
\begin{array}{l}
\displaystyle
\partial_tV+
\tau_b^y+\alpha\bigg(\int_{-H}^{\zeta}
\left(\partial_x(uv)+\partial_yv^2\right))\,dz \\
\\
\displaystyle
\quad-\tau^y_s
-\int_{-H}^{\zeta}\big(
\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle
\quad +fU
-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_y b\,dz'\,dz \bigg)
=
-gD\partial_y\zeta.
\end{array}
\end{equation}

Here, $\tau_b^x$ and $\tau_b^y$ are bottom stresses. Their calculation
is discussed in section \ref{sec-bottom-friction-3d}.
As a first preparation for the mode splitting,
these integrals of the momentum equations can be formally rewritten as

\begin{equation}\label{UMOM}
\begin{array}{l}
\displaystyle
\partial_tU+\frac{R}{D^2}U\sqrt{U^2+V^2}+S^x_{F}
+\alpha\bigg(\partial_x\left(\frac{U^2}{D}\right)
+\partial_y\left(\frac{UV}{D}\right)\\
\\
\displaystyle
\quad-\tau^x_s
-\partial_x\left(2A_h^MD\partial_x\left(\frac{U}{D}\right)\right)
-\partial_y\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right)
\\ \\
\displaystyle
\quad-fV
+S^x_{A}-S^x_D+S^x_B\bigg)
=
-gD\partial_x\zeta
\end{array}
\end{equation}

and

\begin{equation}\label{VMOM}
\begin{array}{l}
\displaystyle
\partial_tV+\frac{R}{D^2}V\sqrt{U^2+V^2}+S^y_{F}
+\alpha\bigg(\partial_x\frac{UV}{D}+\partial_y\frac{V^2}{D}\\
\\
\displaystyle
\quad-\tau^y_s
-\partial_x\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right)
-\partial_y\left(2A_h^MD\partial_y\left(\frac{V}{D}\right)\right)
\\ \\
\displaystyle
\quad+fU
+S^y_{A}-S^y_D+S^y_{B}\bigg)
=
-gD\partial_y\zeta
\end{array}
\end{equation}

with the so-called slow terms
for bottom friction

\begin{equation}\label{Slowfirst}
S^x_{F}=\tau^x_b-\frac{R}{D^2}U\sqrt{U^2+V^2},
\end{equation}

\begin{equation}\label{Slowsecond}
S^y_{F}=\tau^y_b-\frac{R}{D^2}V\sqrt{U^2+V^2},
\end{equation}

horizontal advection

\begin{equation}\label{SxA}
S^x_{A}=\int_{-H}^{\zeta}\left(\partial_x u^2+\partial_y(uv)\right)\,dz-
\partial_x\left(\frac{U^2}{D}\right)-\partial_y\left(\frac{UV}{D}\right),
\end{equation}

\begin{equation}\label{SyA} 
S^y_{A}=\int_{-H}^{\zeta}\left(\partial_x (uv)+\partial_yv^2\right)\,dz-
\partial_x\left(\frac{UV}{D}\right)-\partial_y\left(frac{V^2}{D}\right),
\end{equation}

horizontal diffusion

\begin{equation}\label{SxD}
\begin{array}{l}
\displaystyle
S^x_D=
\int_{-H}^{\zeta}\big(
\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle\qquad 
-\partial_x\left(2A_h^MD\partial_x\left(\frac{U}{D}\right)\right)
-\partial_y\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right), 
\end{array}
\end{equation}

\begin{equation}\label{SyD}
\begin{array}{l}
\displaystyle
S^y_D=
\int_{-H}^{\zeta}\big(
\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle\qquad 
-\partial_y\left(2A_h^MD\partial_y\left(\frac{V}{D}\right)\right)
-\partial_x\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right), 
\end{array}
\end{equation}

and internal pressure gradients

\begin{equation}
S^x_B=-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_x b\,dz'\,dz 
\end{equation}

and

\begin{equation}\label{Slowlast}
S^y_B=-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_y b\,dz'\,dz. 
\end{equation}

The drag coefficient $R$ for the external mode is calculated
as (this logarithmic dependence of the bottom drag from
the water depth and the bottom roughness parameter $z_b^0$ is discussed
in detail by \cite{BURCHARDea02}):


\begin{equation}\label{bottom_vert}
R = \left(\frac{\kappa}
{\ln\left(\frac{\frac{D}{2}+z^b_0}{z^b_0}\right)}\right)^2.
\end{equation}

It should be noted that for numerical reasons, an additional explicit
damping has been implemented into GETM. This method is based on
diffusion of horizontal transports and is described in section

