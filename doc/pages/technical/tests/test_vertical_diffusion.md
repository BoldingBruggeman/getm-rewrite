title: test_vertical_diffusion
Author: Karste Bolding and Hans Burchard

### Vertical diffusion test program

The test program is implemented in [[test_vertical_diffusion.F90]] and in addition to the advection from the [[getm_operators]] implemented as a submodule in [[diffusion.F90]].

#### Test configuration


Using the Gaussian bell function

\begin{equation}
c(z,0) = c_{\max}\exp\left(-\frac{(z-z_0)^2}{a^2}  \right),
\end{equation}

with the height \(a \ll H\) and the vertical position \(z=z_0\) as initial condition, the vertical diffusion equation with diffusivity \(K\) is fulfilled by the following analytical solution:
\begin{equation}
c(z,t)= c_{\max}\sqrt{\frac{a^2}{4K t + a^2}}
\exp\left(-\frac{(z-z_0)^2}{4K t + a^2}  \right).
\end{equation}

