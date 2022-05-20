In GETM, various methods are provided for the calculation of the
internal pressure gradients terms in $x$- and $y$-direction.
These terms which appear as layer-integrated terms in the
equations for the layer-integrated momentum are for the
eastward momentum $p_k$ (see equation (\ref{uEqvi})):

\begin{equation}
h_k\left(\frac12h_N(\partial^*_xb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_xb)_j
\right)
\end{equation}

and for the northward layer-integrated momentum $q_k$
(see equation (\ref{vEqvi})):

\begin{equation}
h_k\left(\frac12h_N(\partial^*_yb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_yb)_j
\right)
\end{equation}

The major problem is how to calculate the horizontal (with respect
to isogeopotentials) buoyancy gradients $\partial^*_xb$ and $\partial^*_yb$,
which need to be defined at the interfaces positioned vertically
between two velocity points.

The methods for calculating the internal pressure gradient included in
GETM are currently:

\begin{enumerate}
\item Method by \cite{MELLORea94}, see routine {\tt ip\_blumberg\_mellor}
\item Modified \cite{MELLORea94} method, exact for linear density profiles
with $z$-dependence only, see routine {\tt ip\_blumberg\_mellor\_lin}
\item Calculation by mean of linear interpolation to $z$-levels, see routine
{\tt ip\_z\_interpol}
\item Method by \cite{SONG98}, see routine {\tt ip\_song\_wright}
\item Method by \cite{CHUea03}, see routine {\tt ip\_chu\_fan}
\item Method by \cite{SHCHEPETKINea03}, see routine {\tt
ip\_shchepetkin\_mcwilliams}
\item Method by \cite{STELLINGea94}, see routine {\tt
ip\_stelling\_vankester.F90}
\end{enumerate}

It is possible, by setting the compiler option {\tt SUBSTR\_INI\_PRESS},
to substract the initial pressure gradient from all pressure
gradients. This is only advisable for strong stratification
without any initial internal pressure gradients. In this case
any non-zero values of the resulting numerical initial pressure gradient
are due to discretisation errors.

