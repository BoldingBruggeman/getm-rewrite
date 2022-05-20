Here, the vertically integrated $V$-momentum equation (\ref{VMOM}) given
on page \pageref{VMOM} including a
number of slow terms is calculated. One slight modification is that
for better stability of drying and flooding processes the slow friction
term $S^y_F$ is now also multiplied with the parameter $\alpha$ defined
in eq.\ (\ref{alpha}).

Furthermore, the horizontal pressure gradient $\partial^*_y\zeta$ is
modified in order to
support drying and flooding, see figure \ref{figpressgrad} on page
\pageref{figpressgrad} and the explanations in section \ref{Section_dry}.
$\partial^*_y\zeta$ is now also considering the atmospheric pressure
gradient at sea surface height.

For numerical stability reasons, the $V$-momentum equation is here
discretised in time such that the
bed friction is treated explicitely:

\begin{equation}\label{Vmom_discrete}
\displaystyle
V^{n+1}=\frac{V^n-\Delta t_m(gD\partial^*_y\zeta
+\alpha(-\frac{\tau_y^s}{\rho_0}+fU^n+V_{Ex}+S_A^y-S_D^y+S_B^y+S_F^y))}
{1+\Delta t_m\frac{R}{D^2}\sqrt{\left(U^n\right)^2+\left(V^n\right)^2}},
\end{equation}

with $V_{Ex}$ combining advection and diffusion of $V$, see routines
{\tt uv\_advect} (section \ref{sec-uv-advect} on page
\pageref{sec-uv-advect}) and {\tt uv\_diffusion}
(section \ref{sec-uv-diffusion} on page
\pageref{sec-uv-diffusion}). The slow terms
are calculated in the routine {\tt slow\_terms} documented in section
\ref{sec-slow-terms} on page \pageref{sec-slow-terms}.
In (\ref{Vmom_discrete}), $V^{n+1}$ denotes the transport on the
new and $U^n$ and $V^n$ the transports on the old time level.

The Coriolis term $fV$ for the subsequent $U$-momentum is also calculated
here, by directly interpolating the $U$-transports to the U-points
or by a method suggested by \cite{ESPELIDea00} which takes the
varying water depths into account.

Some provisions for proper behaviour of the $V$-transports when
GETM runs as slice model are made as well, see section
\ref{Section_GETM_Slice} on page \pageref{Section_GETM_Slice}.

