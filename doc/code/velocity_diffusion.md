title: Velocity diffusion
Author: Karsten Bolding, Hans Burchard and Peter Holterman

Here, the diffusion terms for the vertically integrated transports are
calculated by means of central differences, following the finite volume
approach. They are added to the
advection terms into the terms {\tt UEx} and {\tt VEx} for the
\(U\)- and the \(V\)-equation, respectively. The physical diffusion with the
given eddy viscosity coefficient \(A_h^M\) is based on velocity gradients,
whereas an additional numerical damping of the barotropic mode is based
on gradients of the transports with the damping coefficient \(A_h^N\),
see the example given as equations (\ref{smooth_example_1}) and
(\ref{smooth_example_2}).

First diffusion term in (\ref{UMOM}):
\begin{equation}
\left(mn\,\partial_{\cal X}\left(2A_h^MD\partial_{\cal X}\left(\frac{U}{D}\right)
+ A_h^N\partial_{\cal X}U\right)\right)_{i,j}\approx
\frac{
{\cal F}^{Dxx}_{i+1,j}-{\cal F}^{Dxx}_{i,j}
}{\Delta x^u_{i,j}\Delta y^u_{i,j}}
\end{equation}

with diffusive fluxes

\begin{equation}
{\cal F}^{Dxx}_{i,j}=\left(2 A_h^MD_{i,j}\left(\frac{U_{i,j}}{D^u_{i,j}}
-\frac{U_{i-1,j}}{D^u_{i-1,j}}\right)+A_h^N\left(U_{i,j}
-U_{i-1,j}\right)\right)
\frac{\Delta y^c_{i,j}}{\Delta x^c_{i,j}}.
\end{equation}

Second diffusion term in (\ref{UMOM}):
\begin{equation}
\left(mn\,\partial_{\cal Y}\left(A_h^MD\left(\partial_{\cal Y}\left(\frac{U}{D}\right)+\partial_{\cal X}\left(\frac{V}{D}\right)\right)
+ A_h^N\partial_{\cal Y}U\right)\right)_{i,j}\approx
\frac{
{\cal F}^{Dxy}_{i,j}-{\cal F}^{Dxy}_{i,j-1}
}{\Delta x^x_{i,j}\Delta y^x_{i,j}}
\end{equation}

with diffusive fluxes

\begin{equation}
\begin{array}{rcl}
\displaystyle
{\cal F}^{Dxy}_{i,j}&=&
\displaystyle
A_h^M\frac12\left(D^u_{i,j}+D^u_{i,j+1}\right)
\Delta x^x_{i,j}\left(\left(\frac{U_{i,j+1}}{D^u_{i,j+1}}
-\frac{U_{i,j}}{D^u_{i,j}}\right)\frac{1}{\Delta y^x_{i,j}}
+\left(\frac{V_{i+1,j}}{D^v_{i+1,j}}
-\frac{V_{i,j}}{D^v_{i,j}}\right)\frac{1}{\Delta x^x_{i,j}}\right) \\ \\
&&
\displaystyle
+A_h^N\left(U_{i,j+1} -U_{i,j}\right)\frac{\Delta x^x_{i,j}}
{\Delta y^x_{i,j}}.
\end{array}
\end{equation}

First diffusion term in (\ref{VMOM}):
\begin{equation}
\left(mn\,\partial_{\cal X}\left(A_h^MD\left(\partial_{\cal Y}\left(\frac{U}{D}\right)+\partial_{\cal X}\left(\frac{V}{D}\right)\right)
+ A_h^N\partial_{\cal X}V\right)\right)_{i,j}\approx
\frac{
{\cal F}^{Dyx}_{i,j}-{\cal F}^{Dyx}_{i-1,j}
}{\Delta x^x_{i,j}\Delta y^x_{i,j}}
\end{equation}

with diffusive fluxes

\begin{equation}
\begin{array}{rcl}
\displaystyle
{\cal F}^{Dyx}_{i,j}&=&
\displaystyle
A_h^M\frac12\left(D^v_{i,j}+D^v_{i+1,j}\right)
\Delta y^x_{i,j}\left(\left(\frac{U_{i,j+1}}{D^u_{i,j+1}}
-\frac{U_{i,j}}{D^u_{i,j}}\right)\frac{1}{\Delta y^x_{i,j}}
+\left(\frac{V_{i+1,j}}{D^v_{i+1,j}}
-\frac{V_{i,j}}{D^v_{i,j}}\right)\frac{1}{\Delta x^x_{i,j}}\right) \\ \\
&&
\displaystyle
+A_h^N\left(V_{i+1,j} -V_{i,j}\right)\frac{\Delta y^x_{i,j}}
{\Delta x^x_{i,j}}.
\end{array}
\end{equation}

Second diffusion term in (\ref{VMOM}):
\begin{equation}
\left(mn\,\partial_{\cal Y}\left(2A_h^MD\partial_{\cal Y}\left(\frac{V}{D}\right)
+ A_h^N\partial_{\cal Y}V\right)\right)_{i,j}\approx
\frac{
{\cal F}^{Dyy}_{i,j+1}-{\cal F}^{Dyy}_{i,j}
}{\Delta x^v_{i,j}\Delta y^v_{i,j}}
\end{equation}

with diffusive fluxes

\begin{equation}
{\cal F}^{Dyy}_{i,j}=\left(2 A_h^MD_{i,j}\left(\frac{V_{i,j}}{D^v_{i,j}}
-\frac{V_{i,j-1}}{D^v_{i,j-1}}\right)+A_h^N\left(V_{i,j}
-V_{i,j-1}\right)\right)
\frac{\Delta x^c_{i,j}}{\Delta y^c_{i,j}}.
\end{equation}

The role of the additional diffusion of \(U\) and \(V\) with the
diffusion coefficient \(A_h^N\) is best demonstrated by means of a
simplified set of vertically integrated equations:

\begin{equation}\label{smooth_example_1}
\begin{array}{l}
\displaystyle
\partial_t \eta = - \partial_x U - \partial_y V \\ \\
\displaystyle
\partial_t U = -gD\partial_x\eta
+ A_h^N \left(\partial_{xx} U + \partial_{yy} U\right) \\ \\
\displaystyle
\partial_t V = -gD\partial_y\eta
+ A_h^N \left(\partial_{xx} V + \partial_{yy} V\right), \\ \\
\end{array}
\end{equation}

which can be transformed into an equation for \(\partial_t\eta\) by
derivation of the \(\eta\)-equation with respect to \(t\),
of the \(U\)-equation with respect to \(x\)
and the \(V\)-equation with respect to \(y\) and subsequent elimination of
\(U\) and \(V\):

\begin{equation}\label{smooth_example_2}
\partial_t \left(\partial_t\eta\right) = gD \left(\partial_{xx}\eta
+ \partial_{yy}\eta\right)
+  A_h^N \left(\partial_{xx}\left(\partial_t\eta\right)
+\partial_{yy}\left(\partial_t\eta\right)\right),
\end{equation}

which can be interpreted as a wave equation with a damping on
\(\partial_t \eta\). This introduces an explicit damping of
free surface elevation oscillations in a momentum-conservative
manner. Hydrodynamic models with implicit treatment of the barotropic mode
do not need to apply this method due to the implicit damping of those models,
see e.g.\ \cite{BACKHAUS85}. The implementation of this explicit damping
described here has been suggested by Jean-Marie Beckers, Li\'ege
(Belgium).

When working with the option {\tt SLICE\_MODEL}, the calculation of
all gradients in \(y\)\)-direction is suppressed.

