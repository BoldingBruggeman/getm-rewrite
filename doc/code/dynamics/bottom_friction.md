### Bottom friction

As already mentioned earlier in section \ref{SectionDynBounds},
caution is needed when discretising the
bottom boundary conditions for momentum, (\ref{DynBoundbx}).
They are an example
for a physical condition which has to be modified for the numerical
discretisation, since the discrete velocity point nearest to the bottom
is half a grid box away from the point where the
\index{boundary conditions}boundary condition
is defined. Furthermore, due to the logarithmic \index{law of the wall}
law, high velocity
gradients are typical near the bed. Simply setting the discrete bottom
velocity to zero, would therefore lead to large discretisation errors.
Instead, a flux condition using bottom stresses is derived from the
law of the wall.

For the determination of the normalised bottom stresses

\begin{equation}\label{tauxb}
\frac{\tau^x_b}{\rho_0}=u_*^{bx}u_*^b,
\end{equation}
\begin{equation}\label{tauyb}
\frac{\tau^y_b}{\rho_0}=u_*^{by}u_*^b
\end{equation}

with the friction velocities
$u_*^b=\sqrt{\tau_b/\rho_0}$\label{pagetaub} with
$\tau_b=\sqrt{(\tau^x_{b})^2+(\tau^y_{b})^2}$,
assumptions about the structure of
velocity inside the discrete bottom layer have to be made.
We use here the logarithmic profile

\begin{equation}\label{log_prof}
\frac{u(z')}{u_*}
=\frac{1}{\kappa}\mbox{ln}\left(\frac{z'+z_0^b}{z_0^b}\right),
\end{equation}

with the bottom roughness length $z_0^b$, the von K\'arm\'an constant $\kappa=0.4$
and the distance from the bed, $z'$.
Therefore, estimates for the velocities in the centre of the bottom
layer can be achieved by:

\begin{equation}\label{ulogdis}
u_b = \frac{u_*^{bx}}{\kappa}\mbox{ln} \left(\frac{0.5h_1+z_0^b}{z_0^b}\right),
\end{equation}

\begin{equation}\label{vlogdis}
v_b = \frac{u_*^{by}}{\kappa}\mbox{ln} \left(\frac{0.5h_1+z_0^b}{z_0^b}\right).
\end{equation}

For $h_1\rightarrow 0$, the original \index{boundary conditions!Dirichlet-type}
Dirichlet-type
no-slip boundary conditions (\ref{DynBoundbx}) are retained.
Another possibility would be to specify the bottom velocities $u_b$ and $v_b$
such that they are equal to the layer-averaged \index{law of the wall}
log-law velocities
(see \cite{BAUMERTea92}).
The calculation of this is however slightly more time consuming
and does not lead to a higher accuracy.


