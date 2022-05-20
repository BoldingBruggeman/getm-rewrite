## Mode splitting

The external system consisting of the surface elevation equation
(\ref{Elevation}) and the transport equations (\ref{UMOM}) and
(\ref{VMOM}) underlies a strict time step constraint if the
discretisation is carried out explicitly:

\begin{equation}\label{CourantSurf}
\Delta t < \left[\frac12 \left(\frac{1}{\Delta x}+\frac{1}{\Delta y}\right)
\sqrt{2gD}\right]^{-1}.
\end{equation}

In contrast to that, the time step of the internal system
is only depending on the Courant number for advection,

\begin{equation}\label{CourantVel}
\Delta t < \min\left\{\frac{\Delta x}{u_{\max}},\frac{\Delta y}{v_{\max}}
\right\},
\end{equation}

which in the case of sub-critical flow is a much weaker
constraint. In order not to punish the whole model with a small
time step resulting from the external system, two different approaches
of mode splitting have been developed in the past.

The first approach, in which the external mode is calculated
implicitly, has been proposed by \cite{MADALAea77}.
This method
is numerically stable (if advection is absent)
for unconditionally long time steps.
The temporal approximation is of second order if semi-implicit
treatment is chosen. In such models, the external and internal mode are
generally calculated with the same time steps
(see e.g.\ \cite{BACKHAUS85}). The introduction of
interactions terms like (\ref{Slowfirst})
- (\ref{Slowlast}) is thus not necessary in such models.

Another approach is to use different time steps for the internal
(macro times steps $\Delta t$)
and the external mode \label{dtm} (micro time steps $\Delta t_m$).
One of the first free surface models which has adopted this method is the
Princeton Ocean Model (POM), see \cite{BLUMBERGea87}.
This method has the disadvantage that interaction
terms are needed for the external mode and that the consistency
between internal and external mode is difficult to obtain.
The advantage of this method is that the free surface elevation is
temporally well resolved which is a major requirement for models
including flooding and drying. That is the reason why this
method is adopted here.

The micro time step $\Delta t_m$ has to be an integer fraction $M$ of the
macro time step $\Delta t$. $\Delta t_m$ is limited by the speed of the
surface waves (\ref{CourantSurf}),
$\Delta t$ is limited by the current speed (\ref{CourantVel}).
The time stepping principle is shown in figure \ref{figtimegrid}.
The vertically integrated transports are averaged over each macro time
step:

\begin{equation}\label{Mdef} 
\bar U_{i,j}^{n+1/2} = \frac{1}{M}\sum_{l=n+0.5/M}^{n+(M-0.5)/M} U^l_{i,j}
\end{equation}
and
\begin{equation}
\bar V_{i,j}^{n+1/2} = \frac{1}{M}\sum_{l=n+0.5/M}^{n+(M-0.5)/M} V^l_{i,j}
\end{equation}
such that
\begin{equation}
\begin{array}{l}
\displaystyle
\frac{\zeta_{i,j}^{n+1}-\zeta_{i,j}^{n}}{\Delta t}=
-\frac{\bar U_{i,j}^{n+1/2}-\bar U_{i-1,j}^{n+1/2}}{\Delta x}
-\frac{\bar V_{i,j}^{n+1/2}-\bar V_{i,j-1}^{n+1/2}}{\Delta y}.
\end{array}
\end{equation}


\begin{figure}
\begin{center}
\includegraphics[width=12cm,bbllx=20,bblly=327,bburx=493,bbury=605]{./figures/gridtime.ps}
\caption{Sketch explaining the organisation of the time stepping.
}\label{figtimegrid}
\end{center}
\end{figure}

### Flooding and drying

The main requirement for drying and flooding is that the
vertically integrated fluxes $U$ and $V$ are controlled such
that at no point a negative water depth occurs.
It is clear that parts of the physics which
play an important role in very shallow water of a
few centimetres depth
like non-hydrostatic effects
are not included in the equations.
However, the model is designed in a way that the control of $U$ and $V$
in very shallow water is mainly motivated by the physics
included in the equations
rather than by defining complex drying and flooding algorithms.
It is assumed that the major process in this
situation is a balance between pressure gradient and bottom friction.
Therefore, in the case of very shallow water, all other terms are multiplied
with the non-dimensional factor $\alpha$ which approaches zero
when a minimum water depth is reached.

By using formulation (\ref{bottom_vert})
for calculating the bottom drag coefficient $R$,
it is guaranteed that $R$ is exponentially growing if the water depth
approaches very small values.
This slows the flow down when the water depth in a
velocity point is sinking and also allows for flooding without
further manipulation.

\setlength{\unitlength}{0.5cm}
\begin{figure}
\begin{center}
\fbox{
\begin{picture}(23.0,21.0)(-2.5,-0.5)
\thinlines
\drawline[0](4,8)(4,18)
\drawline[0](10,8)(10,18)
\drawline[0](16,8)(16,18)
\thicklines
\drawline[0](4,15)(10,15)(10,8)(16,8)
\drawline[0](4,15.05)(10.05,15.05)(10.05,8.05)(16,8.05)
\drawline[0](4,14.95)(9.95,14.95)(9.95,7.95)(16,7.95)
\thinlines
\drawline[0](4,16.3)(10,16.3)
\drawline[0](10,11)(16,11)
\dashline[+50]{0.2}(4,16.6)(16,16.6)
\dashline[+50]{0.2}(3,6)(7,6)
\drawline[0](3,4.5)(7,4.5)
\thicklines
\drawline[0](3,3)(7,3)
\drawline[0](3,2.95)(7,2.95)
\drawline[0](3,3.05)(7,3.05)
\put(8,6){\makebox(0,0)[l]{Virtual sea surface elevation}}
\put(8,4.5){\makebox(0,0)[l]{Actual sea surface elevation}}
\put(8,3){\makebox(0,0)[l]{Bathymetry approximation}}
\thinlines
\drawline[10](4,15)(3,14)
\put(2.8,14){\makebox(0,0)[r]{$-H_{i,j}$}}
\drawline[10](4,16.3)(3,16.3)
\put(2.8,16.3){\makebox(0,0)[r]{$\zeta_{i,j}$}}
\drawline[10](4,16.6)(3,17.6)
\put(2.8,17.6){\makebox(0,0)[r]{$-H_{i,j}+H_{\min}$}}
\drawline[10](16,16.6)(17,16.6)
\put(17.2,16.6){\makebox(0,0)[l]{$\tilde \zeta_{i+1,j}$}}
\drawline[10](16,11)(17,11)
\put(17.2,11){\makebox(0,0)[l]{$\zeta_{i+1,j}$}}
\drawline[10](16,8)(17,8)
\put(17.2,8){\makebox(0,0)[l]{$-H_{i+1,j}$}}
\end{picture}
}
\caption{
Sketch explaining the principle of pressure gradient
minimisation during drying and flooding over sloping bathymetry.
}
\label{figpressgrad}
\end{center}
\end{figure}

In this context, one important question is how to calculated the depth
in the velocity points, $H^u$ and $H^v$, since this determines how shallow the
water in the velocity points may become on sloping beaches.
In ocean models, usually, the depth in the velocity points is calculated as
the mean of depths in adjacent elevation points (T-points):

\begin{equation}\label{Hmean}
H^u_{i,j}=\frac12\left(H_{i,j}+H_{i+1,j}\right),
\qquad
H^v_{i,j}=\frac12\left(H_{i,j}+H_{i,j+1}\right).
\end{equation}

Other models which deal with drying and flooding such as the models of
\cite{DUWE88} and \cite{CASULLIea94} use the minimum of the adjacent depths in
the T-points:

\begin{equation}\label{Hmin}
H^u_{i,j}= \min\{H_{i,j},H_{i+1,j}\},
\qquad
H^v_{i,j}= \min\{H_{i,j},H_{i,j+1}\}.
\end{equation}

This guarantees that all depths in the velocity points around a T-point are not
deeper than the depth in the T-point.
Thus, when the T-point depth is approaching
the minimum depth, then all depths in the velocity
points are also small and the
friction coefficient correspondingly large.

Each of the methods has however drawbacks: When the mean is taken as in
(\ref{Hmean}), the risk of negative water depths is relatively big, and thus
higher values of $D_{\min}$ have to be chosen. When the minimum is taken,
large mud-flats might need unrealistically long times for drying since all
the water volume has to flow through relatively shallow velocity boxes.
Also, velocities in these shallow boxes tend to be relatively high in order to
provide sufficient transports. This might lead to numerical instabilities.

Therefore, GETM has both options, (\ref{Hmean}) and (\ref{Hmin})
and the addition of various other options such
as depth depending weighting of the averaging can easily be added.
These options are controlled by the GETM variable {\tt vel\_depth\_method},
see section \ref{sec-uv-depth} (subroutine {\tt uv\_depths}) documented
on page \pageref{sec-uv-depth}.

If a pressure point is dry (i.e.\ its bathymetry value is higher than
a neighbouring sea surface elevation),
the pressure gradient would be unnaturally high with the consequence
of unwanted flow acceleration. Therefore this pressure gradient will
be manipulated such that (only for the pressure gradient calculation)
a virtual sea surface elevation $\tilde \zeta$ is assumed
(see figure \ref{figpressgrad}). In the situation shown in figure
\ref{figpressgrad}, the left pressure point is dry, and the sea surface
elevation
there is for numerical reasons even slightly below the critical value
$-H_{i,j}+H_{\min}$. In order not to let more water flow out of the left cell,
the pressure gradient between the two boxes shown is calculated with a
manipulated sea surface elevation on the right, $\tilde \zeta_{i+1,j}$.

See also \cite{BURCHARDea03a} for a description of drying and flooding
numerics in GETM.

