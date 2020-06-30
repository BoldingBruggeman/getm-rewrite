title: Velocity diffusion
Author: Karsten Bolding, Hans Burchard and Peter Holterman

Here, the shear frequency squared,
\(M^2=\left(\partial_z u\right)^2+\left(\partial_z v\right)^2\)
and the buoyancy frequency squared, \(N^2=\partial_z b\), with buoyancy \(b\) from
(\ref{bdef}) are calculated.
For both calculations, two alternative methods are coded.
The two straight-forward methods which are explained first, do
both have the disadvantage of generating numerical instabilities.
The straight-forward way for calculating \(M^2\) is as follows:

\begin{equation}\label{ShearSquaredOld}
\begin{array}{l}
\displaystyle
(M^2)_{i,j,k}\approx
\frac12
\Bigg(\left(\frac{u_{i,j,k+1}-u_{i,j,k}}
{\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}\right)^2
+
\left(\frac{u_{i-1,j,k+1}-u_{i-1,j,k}}
{\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}\right)^2
\\ \\ \qquad\qquad\quad
\displaystyle
+
\left(\frac{v_{i,j,k+1}-v_{i,j,k}}
{\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}\right)^2
+
\left(\frac{v_{i,j-1,k+1}-v_{i,j-1,k}}
{\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}\right)^2
\Bigg)
\end{array}
\end{equation}
!
\cite{BURCHARD01c} developed a new scheme, which guarantees that
the mean kinetic energy which is dissipated from the mean flow
equals the shear production of turbulent kinetic energy. Therefore,
this scheme should be numerically more stable than (\ref{ShearSquaredOld}):

\begin{equation}\label{ShearSquaredNew}
\begin{array}{l}
\displaystyle
(M^2)_{i,j,k}\approx
\frac12
\Bigg(\frac{\frac12(\nu_{i,j,k}+\nu_{i+1,j,k})
(u_{i,j,k+1}-u_{i,j,k})^2}{\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}
\\ \\ \qquad\qquad\quad
\displaystyle
+
\frac{\frac12(\nu_{i-1,j,k}+\nu_{i,j,k})
(u_{i-1,j,k+1}-u_{i-1,j,k})^2}{\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}
\\ \\ \qquad\qquad\quad
\displaystyle
+
\frac{\frac12(\nu_{i,j,k}+\nu_{i,j+1,k})
(v_{i,j,k+1}-v_{i,j,k})^2}{\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}
\\ \\ \qquad\qquad\quad
\displaystyle
+
\frac{\frac12(\nu_{i,j-1,k}+\nu_{i,j,k})
(v_{i,j-1,k+1}-v_{i,j-1,k})^2}{\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}
\Bigg)
\\ \\ \qquad\qquad\quad
\displaystyle
\cdot
\left(\frac12\left(h^c_{i,j,k}+h^c_{i,j,k+1}\right)\nu_{i,j,k}\right)^{-1}
\end{array}
\end{equation}
!
The straight-forward discretisation of
\(N^2\) is given by
!
\begin{equation}\label{Nstraight}
\begin{array}{l}
\displaystyle
(N^2)_{i,j,k}\approx
\frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}.
\end{array}
\end{equation}

In some cases, together with the straight-forward discretisation
of the shear squared, (\ref{ShearSquaredOld}), this
did not produce stable numerical results. The reason for this might be that
the velocities involved in the calculation for the shear squared do depend
on the buoyancies in the two
neighbouring T-points such that the straight-forward
method (\ref{Nstraight}) leads to an inconsistency.
However, other experiments with the energy-conserving discretisation
of the shear stress squared, (\ref{ShearSquaredNew})
and the straight-forward discretisation of
\(N^2\), (\ref{Nstraight}),  produced numerically stable results.


