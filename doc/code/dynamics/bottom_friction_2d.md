In this routine the bottom friction for the external (vertically integrated)
mode is calculated. This is done separately for the $U$-equation in the
U-points and for the $V$-equation in the V-points.
The drag coefficient $R$ for the external mode is given in eq.\
(\ref{bottom_vert}) on page \pageref{bottom_vert}.
For {\tt runtype=1} (only vertically integrated calculations), the
bottom roughness length is depending on the bed friction
velocity $u_*^b$ and the molecular viscosity $\nu$:

\begin{equation}\label{Defz0b}
z_0^b = 0.1 \frac{\nu}{u_*^b} + \left(z^b_0\right)_{\min},
\end{equation}

see e.g.\ \cite{KAGAN95}, i.e.\ the given roughness may be increased
by viscous effects.
After this, the drag coefficient is multiplied by the absolute value of the
local velocity, which is alculated by dividing the local transports by the
local water depths and by properly interpolating these velocities
to the U- and V-points. The resulting fields are {\tt ru}, representing
$R\sqrt{u^2+v^2}$ on the U-points and {\tt rv}, representing
this quantity on the V-points.

