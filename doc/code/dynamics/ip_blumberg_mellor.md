Here, the internal part of the pressure
gradient
is discretised according to \cite{MELLORea94}.
The crucial part of this term,
which is $(\partial_x^* b)_k$ (in the case of the $u$-equation),
is discretised
between two vertically adjacent velocity points:

\begin{equation}\label{drhodxdiscr}
\begin{array}{l}
\displaystyle
\frac12(h_{i,j,k}+h_{i,j,k+1})\left(m\,\partial_{\cal X}^*b\right)_{i,j,k} \\ \\
\displaystyle
\approx
\frac12(h^u_{i,j,k}+h^u_{i,j,k+1})
\frac{
\frac12 (b_{i+1,j,k+1}+b_{i+1,j,k})-
\frac12 (b_{i,j,k+1}+b_{i,j,k})}
{\Delta x^u_{i,j}}\\ \\
\displaystyle
-
\frac{z^i_{i+1,j,k}-z^i_{i,j,k}}{\Delta x^u_{i,j}}
\left(\frac12 (b_{i+1,j,k+1}+b_{i,j,k+1})-
\frac12 (b_{i+1,j,k}+b_{i,j,k})\right),
\end{array}
\end{equation}

where $z^i_{i,j,k}$ is the $z$-coordinate of the interface in the T-point
above the grid box with the index $(i,j,k)$.

The discretisation of $(\partial_y^* b)_k$ for the $v$-equation is
done accordingly.

In this routine, as a first step, the interface heights are calculated
in the T-points, in order to allow for the calculation of the
coordinate slopes in the U- and V-points. In a second step, the
expression (\ref{drhodxdiscr}) equivalent formulation for the
$y$-direction are integrated up downwards, beginning from the surface.


