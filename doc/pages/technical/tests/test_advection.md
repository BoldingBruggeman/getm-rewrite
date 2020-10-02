title: test_advection
Author: Karste Bolding and Hans Burchard

### Advection test program

The test program is implemented in [[test_advection.F90]] and in addition to the advection from the [[getm_operators]] implemented as a submodule in [[advection.F90]] (and related files with specific advection scheme algorithms - e.g. [[advection_upstream.F90]] via preprocessing [[advection.F90.template]]) also the [[getm_domain]] and [[getm_output]] modules are used.

[[operators.F90]] (open)

[[advection.F90]] (open)


#### Test configuration
This two-dimensional horizontal test case is simulating a solid-body rotation at frequancy \(\omega\), such that the divergence-free velocity field is defined as

\begin{equation*}
u=-\omega y, \quad v=\omega x. 
\end{equation*}

The square domain has a side length of \(L_x=L_y=100~m\) and is discretised by \(100 \times 100\) square grid boxes. The rotational period \(P\) is prescribed (default \(P=600~s\)) from which \(\omega=2\pi/P\) is calculated. Furthermore, the maximum Courant number, \(\mu_{\max}\\) is prescribed (default \(\mu_{\max}=1\)) as well as the number of revolutions, \(M\) (default \(M=5\)).

Since the maximum velocity is \(u_{\max}=\frac12 \omega L_x\), the CFL stability criterium requires
\begin{equation}
\Delta t_{\max}=\mu_{\max} \frac{\Delta x}{u_{\max}}.
\end{equation}

To divide each revolution into an integer number of time steps, the total number of time steps is calculated as 
\begin{equation}
N=M \left\lceil \frac{P}{\Delta t_{\max}} \right\rceil,
\end{equation}

where \(\lceil x\rceil\) is the ceiling function to any non-negative real number \(x\). 

The time step is then calculated as
\begin{equation}
\Delta t = \frac{MP}{N}.
\end{equation}

On the square domain, a scalar field is distributed with an ambient value of \(f=1\) and two solid bodies, the one on the left side being a circular Gaussian bump with a maximum value of \(f=5\) and a width of \(a=10~m\), and the one on the right side being a square body with \(f=5\) of side length \(a=20~m\) with an incision down to the ambient value at one side of \(b=13~m\) length and \(c=6~m\) width.

These two shapes are rotated \(M\) times by advecting them with a choice of the following advection schemes:

  1. HSIMT
  2. MUSCL
  3. P2\_PDM
  4. SPLMAX13
  5. SUPERBEE
  6. UPSTREAM

#### Test results
