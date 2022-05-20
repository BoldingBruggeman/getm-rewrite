Based on the assumption that the velocity distribution in the bottom
layer is logarithmic,
the product of the drag coefficient with the
absolute value of the current speed in the bottom layer,

\begin{equation}
r \sqrt{u_b^2+v_b^2}
\end{equation}

with the velocity components of the bottom layer, $u_b$ and $v_b$,
and the drag coefficient

\begin{equation}\label{r}
r = \left(\frac{\kappa}{\ln \left(\frac{0.5h_1+z_0^b}{z_0^b}\right)}
\right)^2,
\end{equation}

is calculated and
provided as output parameters {\tt rru} (for U-points) and
{\tt rrv} (for V-points). The layer height $h_1$ in (\ref{r}) is set to
the thickness of the bottom layer in the respective U- or V-point.

There are some experimental options for the interested user included
here. It is possible to change the interpolation of $u$ to V-points
and of $v$ to U-points from velocity-based interpolation (as done
presently) to transport-based averaging (commented out). Furthermore,
the user may activate some outcommented lines which allow the
consideration of flow-depending bottom roughness length $z_0^b$
according to (\ref{Defz0b}), see page \pageref{Defz0b}.

For a derivation of (\ref{r}), see section \ref{SectionBedFric} on
page \pageref{SectionBedFric}.

