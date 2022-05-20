Here, the pressure gradient is calculated according to the
method and the algorithm suggested by Shchepetkin and McWilliams,
2003. This method uses a nonconservative  Density-Jacobian scheme,
based on  cubic polynomial fits for the bouyancy "buoy" and "zz",
the vertical position of rho-points, as functions of its respective
array indices. The cubic polynomials are monotonized by using harmonic
mean instead of linear averages to interpolate slopes.
Exact anti-symmetry of the density Jacobian
\begin{equation}
       J(rho,zz)=-J(zz,rho)
\end{equation}
is retained for the density/bouyancy Jacobian in the pressure
gradient formulation in x-direction for a non aligned vertical
coordinate $\sigma$, the atmospheric pressure $p_0$ and the sea
surface elevation $\eta$:
\begin{equation}
\label{eq: shchepetkin pgf}
- {1 \over \rho_0} \partial_x p = \underbrace{-{1 \over \rho_0} \partial_x p_0 - g \partial_x\eta}_{uu\_momentum}
                         + \underbrace{buoy(\eta) \partial_x\eta + \int_z^\eta J(buoy,zz)\mbox{d}\sigma}_{idpdx}
\end{equation}
Details about the calculation of the integral over the Jacobian in
(\ref{eq: shchepetkin pgf}) can be found in Shchepetkin and McWilliams,
2003.

If parameter OneFifth (below) is set to zero, the scheme should
become identical to standard Jacobian.
