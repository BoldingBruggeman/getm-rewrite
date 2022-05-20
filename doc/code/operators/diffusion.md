### General vertical diffusion equation equation

In the following, the discretisation of a simple diffusion equation

\begin{equation}\label{X}
\partial_t X - \partial_z(\nu \partial_z X)=Q
\end{equation}

will be shown for a general physical quantity $X$ which could be momentum $u$ or $v$
or a tracer or turbulent quantity. $Q$ denotes all source terms which will be
discretised on the "old" time level.
After layer integration, (\ref{X}) is of the form

\begin{equation}
\partial_t (h_kX_k) - (\nu \partial_z X)\bigg|_k +(\nu \partial_z X)\bigg|_{k-1} =h_kQ_k
\end{equation}

where $X_k$ and $Q_k$ denote layer averages of $X$ and $Q$, respectively\footnote{
With time-varying depth some extra terms would appear due to the Leibniz rule
(\ref{Leibniz}). These would however be compensated by advective terms which are here
contained in the source term $Q$ for simplicity.}.

with the Neumann-type boundary conditions

 \begin{equation}\label{Fs}
 \nu \partial_z X = F_s
 \qquad \mbox{for } z=\zeta,\qquad
 \end{equation}

 and

 \begin{equation}\label{Fb} 
 \nu \partial_z X = - F_b
 \qquad \mbox{for } z=-H,\qquad
 \end{equation}

the semi-implicit discretisation is for a mean flow quantity
such as momentum or tracers of the following form:

\begin{equation}\label{sigmafirst}
\displaystyle
\frac{h^{n+1}_{N}X^{n+1}_{N}-h^n_{N}X^n_{N}}{\Delta t}
-F_s
+\nu^n_{N-1}\frac{X^{n+\sigma}_{N}-X^{n+\sigma}_{N-1}}{\frac12(h^{n+\sigma}_{N}
+h^{n+\sigma}_{N-1})}
=h^n_NQ^n_N,
\end{equation}

\begin{equation}\label{Xdiscrete}
\begin{array}{l}
\displaystyle
\frac{h^{n+1}_kX^{n+1}_k-h^n_kX^n_k}{\Delta t}
-\nu^n_k\frac{X^{n+\sigma}_{k+1}-X^{n+\sigma}_{k}}{\frac12(h^{n+\sigma}_{k+1}
+h^{n+\sigma}_k)}
\\ \\ \qquad
\displaystyle
+\nu^n_{k-1}\frac{X^{n+\sigma}_{k}-X^{n+\sigma}_{k-1}}{\frac12(h^{n+\sigma}_k
+h^{n+\sigma}_{k-1})}
=h^n_kQ^n_k \qquad \mbox{for } 1<k<N, \qquad
\end{array}
\end{equation}

\begin{equation}
\displaystyle
\frac{h^{n+1}_1X^{n+1}_1-h^n_1X^n_1}{\Delta t}
-\nu^n_1\frac{X^{n+\sigma}_{2}-X^{n+\sigma}_{1}}{\frac12(h^{n+\sigma}_{2}+h^{n+\sigma}_1)}
-F_b
=h^n_1Q^n_1,
\end{equation}

with

\begin{equation}
X_k^{n+\sigma}=\sigma X_k^{n+1}+(1-\sigma)X_k^n
\end{equation}

and

\begin{equation}
h_k^{n+\sigma}=\sigma h_k^{n+1}+(1-\sigma)h_k^n.
\end{equation}
  
Upper indices denote time levels ($n$ for "old" and $n+1$ for "new")
and lower indices denote the vertical discrete location. 
Horizontal indices $(i,j)$ are omitted for simplicity.
Thus, for $\sigma=0$ a fully explicit, for $\sigma=1$ a fully implicit
and for $\sigma=0.5$ the 
Crank-Nicholson second-order in time scheme
is obtained.


This semi-implicit differencing
leads for each transport equation to a system of linear equations
with the following tri-diagonal matrix:
\index{tri-diagonal matrix}

\begin{equation}
\begin{array}{l}
\displaystyle
-X^{n+1}_{N-1} \frac{\sigma\Delta t\, \nu^n_{N-1}}
{\frac12(h^{n+\sigma}_{N}+h^{n+\sigma}_{N-1})}\\ \\
\displaystyle
+X^{n+1}_{N} \left(h^{n+1}_{N}
+\frac{\sigma\Delta t\, \nu^n_{N-1}}{\frac12(h^{n+\sigma}_{N}+h^{n+\sigma}_{{N}-1})}
\right)
= \\ \\
\displaystyle
\qquad X^{n}_{{N}-1} \frac{(1-\sigma)\Delta t\, \nu^n_{{N}-1}}
{\frac12(h^{n+\sigma}_{N}+h^{n+\sigma}_{{N}-1})}\\ \\
\displaystyle
\qquad +X^{n}_{N} \left(h^{n}_{N}
-\frac{(1-\sigma)\Delta t\, \nu^n_{{N}-1}}{\frac12(h^{n+\sigma}_{N}+h^{n+\sigma}_{{N}-1})}
\right)
+\Delta t\,F_s+\Delta t\,h^n_NQ^n_N
\end{array}
\end{equation}

\vspace{1cm}

\begin{equation}
\begin{array}{l}
\displaystyle
-X^{n+1}_{k-1} \frac{\sigma\Delta t\, \nu^n_{k-1}}
{\frac12(h^{n+\sigma}_k+h^{n+\sigma}_{k-1})}\\ \\
\displaystyle
+X^{n+1}_k \left(h_k^{n+1}
+\frac{\sigma\Delta t\, \nu^n_{k-1}}{\frac12(h^{n+\sigma}_k+h^{n+\sigma}_{k-1})}
+\frac{\sigma\Delta t\, \nu^n_{k  }}{\frac12(h^{n+\sigma}_{k+1}+h^{n+\sigma}_k)}
\right)\\ \\
\displaystyle
-X^{n+1}_{k+1} \frac{\sigma\Delta t\, \nu^n_k}
{\frac12(h^{n+\sigma}_{k+1}+h^{n+\sigma}_k)}
= \\ \\
\displaystyle
\qquad X^{n}_{k-1} \frac{(1-\sigma)\Delta t\, \nu^n_{k-1}}
{\frac12(h^{n+\sigma}_k+h^{n+\sigma}_{k-1})}\\ \\
\displaystyle
\qquad +X^{n}_k \left(h_k^n
-\frac{(1-\sigma)\Delta t\, \nu^n_{k-1}}{\frac12(h^{n+\sigma}_k+h^{n+\sigma}_{k-1})}
-\frac{(1-\sigma)\Delta t\, \nu^n_{k  }}{\frac12(h^{n+\sigma}_{k+1}+h^{n+\sigma}_k)}
\right)\\ \\
\displaystyle
\qquad +X^{n}_{k+1} \frac{(1-\sigma)\Delta t\, \nu^n_k}
{\frac12(h^{n+\sigma}_{k+1}+h^{n+\sigma}_k)} +\Delta t\,h^n_kQ^n_k
\qquad \mbox{for } 1<k<{N}, \qquad
\end{array}
\end{equation}

\vspace{1cm}

\begin{equation}
\begin{array}{l}
\displaystyle
X^{n+1}_1 \left(h^{n+1}_1
+\frac{\sigma\Delta t\, \nu^n_{1  }}{\frac12(h^{n+\sigma}_{2}+h^{n+\sigma}_1)}
\right)\\ \\
\displaystyle
-X^{n+1}_{2} \frac{\sigma\Delta t\, \nu^n_1}
{\frac12(h^{n+\sigma}_{2}+h^{n+\sigma}_1)}
=\\ \\
\displaystyle
\qquad X^{n}_1 \left(h^{n}_1
-\frac{(1-\sigma)\Delta t\, \nu^n_{1  }}{\frac12(h^{n+\sigma}_{2}+h^{n+\sigma}_1)}
\right)\\ \\
\displaystyle
\qquad +X^{n}_{2} \frac{(1-\sigma)\Delta t\, \nu^n_1}
{\frac12(h^{n+\sigma}_{2}+h^{n+\sigma}_1)}
+\Delta t\,F_b+\Delta t\,h^n_1Q^n_1,
\end{array}
\end{equation}

which is solved by means of the simplified
Gaussian elimination, see e.g.\ \cite{SAMARSKIJ84}.

It should be noted that the source term $Q$ is often treated
quasi-implicitly,
following a suggestion made by \cite{PATANKAR80}
(see also \cite{DELEERSNIJDERea97}).  
This is done in order to guarantee positive discrete solutions for
physically positive quantities such as concentrations or
turbulent quantities. 
For details of this numerical scheme, which used for
positive quantities in GETM, see also \cite{BURCHARD01b}. 

