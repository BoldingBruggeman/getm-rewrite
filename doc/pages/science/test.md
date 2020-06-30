title: Testing
author: Karsten Bolding

****

General Estuarine Transport Model

![image](./figures/GETM2.ps){width="\linewidth"}

****

Source Code and Test Case\
Documentation

****

Hans Burchard$^{1}$, Karsten Bolding$^2$\
and Lars Umlauf$^1$

****

1: Baltic Sea Research Institute Warnemünde, Germany\
2: Bolding & Bruggeman ApS, Asperup, Denmark\

****

Version pre\_2.6.x (based on v2.5.x code)

****

What’s new
==========

-   [2011-04-11 - 2012-04-01: Development of stable version v2.2]{}

Introduction
============

What is GETM?
-------------

A short history of GETM {#SectionHistory}
-----------------------

The idea for GETM was born in May 1997 in Arcachon, France during a
workshop of the PhaSE project which was sponsored by the European
Community in the framework of the MAST-III programme. It was planned to
set up an idealised numerical model for the Eastern Scheldt, The
Netherlands for simulating the effect of vertical mixing of nutrients on
filter feeder growth rates. A discussion between the first author of
this report, Peter Herman (NIOO, Yerseke, The Netherlands) and Walter
Eifler (JRC Ispra, Italy) had the result that the associated processes
were inherently three-dimensional (in space), and thus, only a
three-dimensional model could give satisfying answers. Now the question
arose, which numerical model to use. An old wadden sea model by
[@BURCHARD95] including a two-equation turbulence model was written in
$z$-coordinates with fixed geopotential layers (which could be added or
removed for rising and sinking sea surface elevation, respectively) had
proven to be too noisy for the applications in mind. Furthermore, the
step-like bottom approximation typical for such models did not seem to
be sufficient. Other Public Domain models did not allow for drying and
flooding of inter-tidal flats, such as the Princeton Ocean Model (POM).
There was thus the need for a new model. Most of the ingredients were
however already there. The first author of this report had already
written a $k$-$\eps$ turbulence model, see [@BURCHARDea95], the
forerunner of GOTM. A two-dimensional code for general vertical
coordinates had been written as well, see [@BURCHARDea97]. And the first
author of this report had already learned a lot about mode splitting
models from Jean-Marie Beckers (University of Liege, Belgium). Back from
Arcachon in Ispra, Italy at the Joint Research Centre of the European
Community, the model was basically written during six weeks, after which
an idealised tidal simulation for the Sylt-Rømø Bight in the wadden sea
area between Germany and Denmark could be successfully simulated, see
[@BURCHARD98]. By that time this model had the little attractive name
[*MUDFLAT*]{} which at least well accounted for the models ability to
dry and flood inter-tidal flats. At the end of the PhaSE project in
1999, the idealised simulation of mussel growth in the Eastern Scheldt
could be finished (not yet published, pers. comm. Francois Lamy and
Peter Herman).

In May 1998 the second author of this report joined the development of
[*MUDFLAT*]{}. He first fully rewrote the model from a one-file
FORTRAN77 code to a modular FORTRAN90/95 code, made the interface to
GOTM (such that the original $k$-$\eps$ model was not used any more),
integrated the netCDF-library into the model, and prepared the
parallelisation of the model. And a new name was created, GETM, General
Estuarine Transport Model. As already in GOTM, the word “General” does
not imply that the model is general, but indicates the motivation to
make it more and more general.

At that time, GETM has actually been applied for simulating currents
inside the Mururoa atoll in the Pacific Ocean, see [@MATHIEUea01].

During the year 2001, GETM was then extended by the authors of this
report to be a fully baroclinic model with transport of active and
passive tracers, calculation of density, internal pressure gradient and
stratification, surface heat and momentum fluxes and so forth. During a
stay of the first author at the Université Catholique de Louvain,
Institut d’Astronomie et de Géophysique George Lemaître, Belgium (we are
grateful to Eric Deleersnijder for this invitation and many discussions)
the high-order advection schemes have been written. During another
invitation to Belgium, this time to the GHER at the Université de Liège,
the first author had the opportunity to discuss numerical details of
GETM with Jean-Marie Beckers, who originally motivated us to use the
mode splitting technique.

The typical challenging application in mind of the authors was always a
simulation of the tidal Elbe, where baroclinicity and drying and
flooding of inter-tidal flats play an important role. Furthermore, the
tidal Elbe is long, narrow and bended, such that the use of Cartesian
coordinates would require an indexing of the horizontal fields, see
e.g. [@DUWE88]. Thus, the use of curvi-linear coordinates which follow
the course of the river has already been considered for a long time.
However, the extensions just listed above, give the model also the
ability to simulate shelf sea processes in fully baroclinic mode, such
that the name General Estuarine Transport Model is already a bit too
restrictive.

The physical equations behind GETM
==================================

Hydrodynamic equations
----------------------

### Three-dimensional momentum equations {#Section_3d_momentum}

For geophysical coastal sea and ocean dynamics, usually the
three-dimensional hydrostatic equations of motion with the Boussinesq
approximation and the eddy viscosity assumption are used ([@BRYAN69],
[@COX84], [@BLUMBERGea87], [@HAIDVOGELea99], [@KANTHAea00b]). In the
flux form, the dynamic equations of motion for the horizontal velocity
components can be written in Cartesian coordinates as:

$$\label{uEq}
\begin{array}{l}
\displaystyle
\partial_t u
+\partial_z(uw)
-\partial_z\left((\nu_t+\nu) \partial_z u\right)
\\ \\ \displaystyle
\qquad+\alpha\bigg(\partial_x(u^2)+\partial_y(uv)
-\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)
\\ \\ \displaystyle
\qquad
-fv
-\int_z^{\zeta}\partial_x b\,dz' \bigg)
=
- g\partial_x \zeta,
\end{array}$$

$$\label{vEq}
\begin{array}{l}
\displaystyle
\partial_t v +\partial_z(vw)
-\partial_z\left((\nu_t+\nu) \partial_z v\right)
\\ \\
\displaystyle
\qquad+\alpha\bigg(\partial_x(vu)+\partial_y(v^2)
-\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)
\\ \\
\displaystyle
\qquad
+fu
-\int_z^{\zeta}\partial_y b\,dz' \bigg)=
- g\partial_y \zeta.
\end{array}$$

The vertical velocity is calculated by means of the incompressibility
condition:

$$\label{Konti}
\partial_x u +\partial_y v +\partial_z w = 0.$$

Here, $u$, $v$ and $w$ are the ensemble averaged velocity components
with respect to the $x$, $y$ and $z$ direction, respectively. The
vertical coordinate $z$ ranges from the bottom $-H(x,y)$ \[Hxy\] to the
surface $\zeta(t,x,y)$ with $t$ denoting time. $\nu_t$ is the vertical
eddy viscosity, $\nu$ the kinematic viscosity, $f$ the Coriolis
parameter, and $g$ is the gravitational acceleration. The horizontal
mixing is parameterised by terms containing the horizontal eddy
viscosity $A_h^M$, see [@BLUMBERGea87]. The buoyancy $b$ is defined as

$$\label{bdef}
b=-g\frac{\rho-\rho_0}{\rho_0}$$

with the density $\rho$ and a reference density $\rho_0$. The last term
on the left hand sides of equations (\[uEq\]) and (\[vEq\]) are the
internal (due to density gradients) and the terms on the right hand
sides are the external (due to surface slopes) pressure gradients. In
the latter, the deviation of surface density from reference density is
neglected (see [@BURCHARDea97]). The derivation of equations (\[uEq\]) -
(\[Konti\]) has been shown in numerous publications, see
e.g. [@PEDLOSKY87], [@HAIDVOGELea99], [@BURCHARD02].

In hydrostatic 3D models, the vertical velocity is calculated by means
of equation (\[Konti\]) velocity equation. Due to this, mass
conservation and free surface elevation can easily be obtained.

Drying and flooding of mud-flats is already incorporated in the physical
equations by multiplying some terms with the non-dimensional number
$\alpha$ which equals unity in regions where a critical water depth
$D_{crit}$ is exceeded and approaches zero when the water depth $D$
tends to a minimum value $D_{min}$:

$$\label{alpha}
\alpha=\min\left\{1,\frac{D-D_{min}}{D_{crit}-D_{min}}\right\}.$$

Thus, $\alpha=1$ for $D\geq D_{crit}$, such that the usual momentum
equation results except for very shallow water, where simplified physics
are considered with a balance between tendency, friction and external
pressure gradient. In a typical wadden sea application, $D_{crit}$ is of
the order of 0.1 m and $D_{\min}$ of the order of 0.02 m (see
[@BURCHARD98], [@BURCHARDea03a]).

### Kinematic boundary conditions and surface elevation equation

At the surface and at the bottom, kinematic boundary conditions result
from the requirement that the particles at the boundaries are moving
along these boundaries:

$$\label{KinBoundsSurf}
w= \partial_t \zeta +  u \partial_x\zeta + v \partial_y \zeta
\qquad \mbox{for } z=\zeta,$$

$$\label{KinBoundsBott}
w= - u \partial_xH - v \partial_y H
\qquad \mbox{for } z=-H.$$

### Dynamic boundary conditions {#SectionDynBounds}

At the bottom boundaries, no-slip conditions are prescribed for the
horizontal velocity components:

$$\label{DynBoundbx}
u=0, \quad v=0.$$

With (\[KinBoundsBott\]), also $w=0$ holds at the bottom. It should be
noted already here, that the bottom boundary condition (\[DynBoundbx\])
is generally not directly used in numerical ocean models, since the
near-bottom values of the horizontal velocity components are not located
at the bed, but half a grid box above it. Instead, a logarithmic
velocity profile is assumed in the bottom layer, leading to a quadratic
friction law, see section \[sec-bottom-friction-3d\].

At the surface, the dynamic boundary conditions read:

$$\label{DynBoundSx}
\begin{array}{l}
(\nu_t+\nu) \partial_z u=\alpha \tau^x_{s},
\\ \\
(\nu_t+\nu) \partial_z v=\alpha \tau^y_{s},
\end{array}$$

The surface stresses (normalised by the reference density) $\tau_s^x$
and $\tau_s^y$ are calculated as functions of wind speed, wind
direction, surface roughness etc. Also here, the drying parameter
$\alpha$ is included in order to provide an easy handling of drying and
flooding.

### Lateral boundary conditions

Let $G$\[G\] denote the lateral boundary of the model domain with the
closed land boundary $G^c$\[Gc\] and the open boundary $G^o$\[Go\] such
that $G^c \cup G^o=G$ and $G^c \cap G^o=\emptyset$. Let further
$\vec u=(u,v)$\[vecu\] denote the horizontal velocity vector and
$\vec u_n=(-v,u)$\[vecun\] its normal vector. At closed boundaries, the
flow must be parallel to the boundary:

$$\vec u_n \cdot \vec\nabla G^c = 0$$

with $\vec\nabla=(\partial_x,\partial_y)$ being the gradient operator.

For an eastern or a western closed boundary with $\vec\nabla G^c=(0,1)$
this has the consequence that $u=0$ and, equivalently, for a southern or
a northern closed boundary with $\vec\nabla G^c=(1,0)$ this has the
consequence that $v=0$.

At open boundaries, the velocity gradients across the boundary vanish:

$$\vec\nabla_n u \cdot \vec\nabla G^o = 0, \quad 
\vec\nabla_n v \cdot \vec\nabla G^o = 0, \quad$$

with $\vec\nabla_n=(-\partial_y,\partial_x)$ being the operator normal
to the gradient operator.

For an eastern or a western open boundary with this has the consequence
that $\partial_x u=\partial_x v =0$ and, equivalently, for a southern or
a northern open boundary this has the consequence that
$\partial_y u=\partial_y v =0$.

At so-called forced open boundaries, the sea surface elevation $\zeta$
is prescribed. At passive open boundaries, it is assumed that the
curvature of the surface elevation normal to the boundary is zero, with
the consequence that the spatial derivatives of the surface slopes
normal to the boundaries vanish.

GETM as slice model {#Section_GETM_Slice}
-------------------

By chosing the compiler option [SLICE\_MODEL]{} it is possible to
operate GETM as a two-dimensional vertical ($xz$-)model under the
assumption that all gradients in $y$-direction vanish. In order to do
so, a bathymetry file with a width of 4 grid points has to be generated,
with the outer ($j=1$, $j=4$) bathymetry values set to land, and the two
inner ones being independent on $j$. The compiler option
[SLICE\_MODEL]{} then sets the transports, velocities, and sea surface
elevations such that they are independent of $y$, i.e. they are forced
to be identical for the same $j$-index. Especially, the $V$-transports
and velocities in the walls ($j=1$, $j=3$) are set to the calculated
value at index $j=2$.

Transformations {#SectionTransformations}
===============

General vertical coordinates {#SectionGeneralCoordinates}
----------------------------

As a preparation of the discretisation, the physical space is vertically
divided into $N$ layers. This is done by introducing internal surfaces
$z_k$, $k=1,\dots,N-1$ which do not intersect, each depending on the
horizontal position $(x,y)$ and time $t$. Let

$$\label{gamma_def}
-H(x,y)=z_0(x,y)<z_1(x,y,t)<\dots<z_{N-1}(x,y,t)<z_N(x,y,t)=\zeta(x,y,t)$$

define the local layer depths $h_k$ with

$$\label{hkdef}
h_k=z_k-z_{k-1}.$$

for $1\leq k\leq N$. For simplicity, the argument $(x,y,t)$ is omitted
in most of the cases.

The most simple layer distribution is given by the so-called $\sigma$
\[sigmacoord\] transformation (see [@PHILLIPS57] for a first application
in meteorology and [@FREEMANea72] for a first application in
hydrodynamics) with

$$\label{sigma}
\sigma_k=\frac{k}{N}-1$$

and

$$\label{sigmaz}
z_k=D\sigma_k$$

for $0\leq k\leq N$.

The $\sigma$-coordinates can also be refined towards the surface and the
bed:

$$\label{formula_Antoine}
\beta_k = \frac{\mbox{tanh}\left( (d_l+d_u)(1+\sigma_k)-d_l\right)
+\mbox{tanh}(d_l)}{\mbox{tanh}(d_l)+\mbox{tanh}(d_u)}-1,
\qquad k=0,\dots,N\qquad$$

such that $z$-levels are obtained as follows:

$$z_k=D\beta_k$$

for $0\leq k\leq N$.

The grid is refined towards the surface for $d_u>0$ and refined towards
the bottom for $d_l>0$. When both, $d_u$ and $d_l$ are larger than zero,
then refinement towards surface and bed is obtained. For $d_u=d_l=0$ the
$\sigma$-transformation (\[sigma\]) with $\beta_k=\sigma_k$ is retained.
Figure \[FigGeneral1\] shows four examples for vertical layer
distributions obtained with the $\sigma$-transformation.

Due to the fact that all layer thicknesses are proportional to the water
depth, the equidistant and also the non-equidistant
$\sigma$-transformations, (\[sigma\]) and (\[formula\_Antoine\]), have
however one striking disadvantage. In order to sufficiently resolve the
mixed layer also in deep water, many layers have to be located near the
surface. The same holds for the bottom boundary layer. This problem of
$\sigma$-coordinates has been discussed by several authors (see
e.g. [@DELEERSNIJDERea92], [@KOK92], [@GERDES93], [@SONGea94],
[@BURCHARDea97]) who suggested methods for generalised vertical
coordinates not resulting in layer thicknesses not proportional to the
water depth.

The generalised vertical coordinate introduced here is a generalisation
of the so-called mixed-layer transformation suggested by
[@BURCHARDea97]. It is a hybrid coordinate which interpolates between
the equidistant and the non-equidistant $\sigma$-transformations given
by (\[sigma\]) and (\[formula\_Antoine\]). The weight for the
interpolation depends on the ratio of a critical water depth
$D_{\gamma}$ (below which equidistant $\sigma$-coordinates are used) and
the actual water depth:

$$\label{HJVTrans}
z_k = D\left(\alpha_{\gamma} \sigma_k + (1-\alpha_{\gamma}) \beta_k\right)$$

with

$$\label{MLDTransform}
\alpha_{\gamma} = \min\left(\frac{
(\beta_k-\beta_{k-1})-\frac{D_{\gamma}}{D}(\sigma_k-\sigma_{k-1})}
{(\beta_k-\beta_{k-1})-(\sigma_k-\sigma_{k-1})},1\right).$$

and $\sigma_k$ from (\[sigma\]) and $\beta_k$ from
(\[formula\_Antoine\]).

For inserting $k=N$ in (\[MLDTransform\]) and $d_l=0$ and $d_u>0$ in
(\[formula\_Antoine\]), the mixed layer transformation of
[@BURCHARDea97] is retained, see the upper two panels in figure
\[FigGeneral2\]. Depending on the values for $D_{\gamma}$ and $d_u$,
some near-surface layer thicknesses will be constant in time and space,
allowing for a good vertical resolution in the surface mixed layer.

The same is obtained for the bottom with the following settings: $k=1$,
$d_l>0$ and $d_u=0$, see the lower two panels in figure \[FigGeneral2\].
This is recommended for reproducing sedimentation dynamics and other
benthic processes. For $d_l=d_u>0$ and $k=1$ or $k=N$ a number of layers
near the surface and near the bottom can be fixed to constant thickness.
Intermediate states are obtained by intermediate settings, see figure
\[FigGeneral3\]. Some pathological settings are also possible, such as
$k=1$, $d_l=1.5$ and $d_u=5$, see figure \[FigGeneral4\].

\[cc\]\[\]\[0.6\][$x$ / nm]{} \[cc\]\[\]\[0.6\][$z$ / m]{}
\[cc\]\[\]\[0.8\][$\sigma$-coordinates, $d_u=0$, $d_l=0$]{} ![
$\sigma$-transformation with four different zooming options. The plots
show the vertical layer distribution for a cross section through the
North Sea from Scarborough in England to Esbjerg in Denmark. The shallow
area at about $x=100$ nm is the Doggerbank.
[]{data-label="FigGeneral1"}](./figures/sigma.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][$\sigma$-coordinates, $d_u=1.5$, $d_l=0$]{} ![
$\sigma$-transformation with four different zooming options. The plots
show the vertical layer distribution for a cross section through the
North Sea from Scarborough in England to Esbjerg in Denmark. The shallow
area at about $x=100$ nm is the Doggerbank.
[]{data-label="FigGeneral1"}](./figures/beta10.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][$\sigma$-coordinates, $d_u=0$, $d_l=1.5$]{} ![
$\sigma$-transformation with four different zooming options. The plots
show the vertical layer distribution for a cross section through the
North Sea from Scarborough in England to Esbjerg in Denmark. The shallow
area at about $x=100$ nm is the Doggerbank.
[]{data-label="FigGeneral1"}](./figures/beta01.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][$\sigma$-coordinates, $d_u=1.5$, $d_l=1.5$]{} ![
$\sigma$-transformation with four different zooming options. The plots
show the vertical layer distribution for a cross section through the
North Sea from Scarborough in England to Esbjerg in Denmark. The shallow
area at about $x=100$ nm is the Doggerbank.
[]{data-label="FigGeneral1"}](./figures/beta11.ps "fig:"){width="7cm"}

\[cc\]\[\]\[0.6\][$x$ / nm]{} \[cc\]\[\]\[0.6\][$z$ / m]{}
\[cc\]\[\]\[0.8\][upper $\gamma$-coordinates, $d_u=1.5$, $d_l=0$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in the surface mixed layer (upper two panels)
and with concentration of layers in the bottom mixed layer (lower two
panels). The critical depth $D_{\gamma}$ is here set to 20 m, such that
at all shallower depths the equidistant $\sigma$-transformation is used.
The same underlying bathymetry as in figure \[FigGeneral1\] has been
used.
[]{data-label="FigGeneral2"}](./figures/gammaup1.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][upper $\gamma$-coordinates, $d_u=5$, $d_l=0$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in the surface mixed layer (upper two panels)
and with concentration of layers in the bottom mixed layer (lower two
panels). The critical depth $D_{\gamma}$ is here set to 20 m, such that
at all shallower depths the equidistant $\sigma$-transformation is used.
The same underlying bathymetry as in figure \[FigGeneral1\] has been
used.
[]{data-label="FigGeneral2"}](./figures/gammaup5.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][lower $\gamma$-coordinates, $d_u=0$, $d_l=1.5$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in the surface mixed layer (upper two panels)
and with concentration of layers in the bottom mixed layer (lower two
panels). The critical depth $D_{\gamma}$ is here set to 20 m, such that
at all shallower depths the equidistant $\sigma$-transformation is used.
The same underlying bathymetry as in figure \[FigGeneral1\] has been
used.
[]{data-label="FigGeneral2"}](./figures/gammalow1.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][lower $\gamma$-coordinates, $d_u=0$, $d_l=5$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in the surface mixed layer (upper two panels)
and with concentration of layers in the bottom mixed layer (lower two
panels). The critical depth $D_{\gamma}$ is here set to 20 m, such that
at all shallower depths the equidistant $\sigma$-transformation is used.
The same underlying bathymetry as in figure \[FigGeneral1\] has been
used.
[]{data-label="FigGeneral2"}](./figures/gammalow5.ps "fig:"){width="7cm"}

\[cc\]\[\]\[0.6\][$x$ / nm]{} \[cc\]\[\]\[0.6\][$z$ / m]{}
\[cc\]\[\]\[0.8\][upper $\gamma$-coordinates, $d_u=5$, $d_l=1.5$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in both, the surface mixed layer and the bottom
mixed layer. Four different realisations are shown. The critical depth
$D_{\gamma}$ is here set to 20 m, such that at all shallower depths the
equidistant $\sigma$-transformation is used. The same underlying
bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral3"}](./figures/gammauplow.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][lower $\gamma$-coordinates, $d_u=1.5$, $d_l=5$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in both, the surface mixed layer and the bottom
mixed layer. Four different realisations are shown. The critical depth
$D_{\gamma}$ is here set to 20 m, such that at all shallower depths the
equidistant $\sigma$-transformation is used. The same underlying
bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral3"}](./figures/gammalowup.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][symmetric $\gamma$-coordinates, $d_u=1.5$,
$d_l=1.5$]{} ![ Boundary layer transformation (or $\gamma$
transformation) with concentration of layers in both, the surface mixed
layer and the bottom mixed layer. Four different realisations are shown.
The critical depth $D_{\gamma}$ is here set to 20 m, such that at all
shallower depths the equidistant $\sigma$-transformation is used. The
same underlying bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral3"}](./figures/gammauplow1.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][symmetric $\gamma$-coordinates, $d_u=5$, $d_l=5$]{} ![
Boundary layer transformation (or $\gamma$ transformation) with
concentration of layers in both, the surface mixed layer and the bottom
mixed layer. Four different realisations are shown. The critical depth
$D_{\gamma}$ is here set to 20 m, such that at all shallower depths the
equidistant $\sigma$-transformation is used. The same underlying
bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral3"}](./figures/gammauplow5.ps "fig:"){width="7cm"}

\[cc\]\[\]\[0.6\][$x$ / nm]{} \[cc\]\[\]\[0.6\][$z$ / m]{}
\[cc\]\[\]\[0.8\][upper $\gamma$-coordinates, $d_u=1.5$, $d_l=5$]{} ![
Two pathological examples for the boundary layer transformation. The
critical depth $D_{\gamma}$ is here set to 20 m, such that at all
shallower depths the equidistant $\sigma$-transformation is used. The
same underlying bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral4"}](./figures/gammapathoup.ps "fig:"){width="7cm"}
\[cc\]\[\]\[0.8\][lower $\gamma$-coordinates, $d_u=5$, $d_l=1.5$]{} ![
Two pathological examples for the boundary layer transformation. The
critical depth $D_{\gamma}$ is here set to 20 m, such that at all
shallower depths the equidistant $\sigma$-transformation is used. The
same underlying bathymetry as in figure \[FigGeneral1\] has been used.
[]{data-label="FigGeneral4"}](./figures/gammapatholow.ps "fig:"){width="7cm"}

The strong potential of the general vertical coordinates concept is the
extendibility towards vertically adaptive grids. Since the layers may be
redistributed after every baroclinic time step, one could adapt the
coordinate distribution to the internal dynamics of the flow. One could
for example concentrate more layers at vertical locations of high
stratification and shear, or force certain layer interfaces towards
certain isopycnals, or approximate Lagrangian vertical coordinates by
minimising the vertical advection through layer interfaces. The
advantages of this concept have recently been demonstrated for
one-dimensional water columns by [@BURCHARDea04]. The three-dimensional
generalisation of this concept of adaptive grids for GETM is currently
under development.

Layer-integrated equations {#SectionLayerIntegrated}
--------------------------

There are two different ways to derive the layer-integrated equations.
[@BURCHARDea97] transform first the equations into general vertical
coordinate form (see [@DELEERSNIJDERea92]) and afterwards integrate the
transformed equations over constant intervals in the transformed space.
[@LANDERea94] integrate the equations in the Cartesian space over
surfaces $z_k$ by considering the Leibniz rule

$$\label{Leibniz}
\int_{z_{k-1}}^{z_k}\partial_x f\,dz
=
\partial_x\int_{z_{k-1}}^{z_k} f\,dz
-f(z_k)\partial_xz_k
+f(z_{k-1})\partial_xz_{k-1}$$

for any function $f$. For the vertical staggering of the layer notation
see figure \[figvertgrid\].

More details about the layer integration are given in [@BURCHARDea97].

With the further definitions of layer integrated transport,

$$\label{pqdef} 
p_k:=\int_{z_{k-1}}^{z_k}u\,dz,\qquad
q_k:=\int_{z_{k-1}}^{z_k}v\,dz,$$

layer mean velocities,

$$\label{ukvkdef} 
u_k:=\frac{p_k}{h_k},\qquad v_k:=\frac{q_k}{h_k},$$

and layer averaged tracer concentrations and buoyancy,

$$\label{ckbkdef} 
c^i_k:=\frac{1}{h_k}\int_{z_{k-1}}^{z_k}c^i\,dz,\qquad
b_k:=\frac{1}{h_k}\int_{z_{k-1}}^{z_k}b\,dz,$$

and the grid related vertical velocity,

$$\label{barwdef} 
\bar w_k:=(w-\partial_tz-u\partial_xz-v\partial_yz)_{z=z_k},$$

the continuity equation (\[Konti\]) has the layer-integrated form:

$$\label{ContiLayerInt}
\partial_t h_k + \partial_x p_k + \partial_y q_k + \bar w_k - \bar w_{k-1}=0.$$

It should be noted that the grid related velocity is located on the
layer interfaces. After this, the layer-integrated momentum equations
read as:

$$\label{uEqvi}
\begin{array}{l}
\partial_t p_k 
+\bar w_k \tilde u_k -\bar w_{k-1} \tilde u_{k-1} 
-\tau^x_k + \tau^x_{k-1} 
\\ \\ \quad
+\alpha\Bigg\{\partial_x(u_kp_k)+\partial_y(v_kp_k)
\\ \\ \displaystyle \quad 
-\partial_x\left(2A_k^Mh_k\partial_xu_k\right)-\partial_y\left(A_k^Mh_k
(\partial_yu_k+\partial_xv_k)\right)
- fq_k 
\\ \\ \quad
\displaystyle
-h_k\left(
\frac12h_N(\partial^*_xb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_xb)_j
\right)\Bigg\}
=
-gh_k\partial_x\zeta,
\end{array}$$

$$\label{vEqvi}
\begin{array}{l}
\partial_t q_k 
+\bar w_k  \tilde v_k -\bar w_{k-1}  \tilde v_{k-1} 
-\tau^y_k + \tau^y_{k-1} 
\\ \\ \quad
+\alpha\Bigg\{\partial_x(u_kq_k)+\partial_y(v_kq_k)
\\ \\ \displaystyle \quad 
-\partial_y\left(2A_k^Mh_k\partial_yv_k\right)-\partial_x\left(A_k^Mh_k
(\partial_yu_k+\partial_xv_k)\right)
+ fp_k 
\\ \\ \quad
\displaystyle
-h_k\left(
\frac12h_N(\partial^*_yb)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_yb)_j
\right)\Bigg\}
=
-gh_k\partial_y\zeta
\end{array}$$

with suitably chosen advective horizontal velocities $\tilde u_k$ and
$\tilde v_k$ (see section \[sec-uv-advect-3d\]) on page , the shear
stresses

$$\label{tauxkdef}
\tau^x_k = \left(\nu_t \partial_z u \right)_k,$$

and

$$\label{tauykdef}
\tau^y_k = \left(\nu_t \partial_z v \right)_k,$$

and the horizontal buoyancy gradients

$$(\partial^*_xb)_k=\frac12(\partial_xb_{k+1}+\partial_x b_k)
-\partial_xz_k\frac{b_{k+1}-b_k}{\frac12(h_{k+1}+h_k)}$$

and

$$(\partial^*_yb)_k=\frac12(\partial_yb_{k+1}+\partial_y b_k)
-\partial_yz_k\frac{b_{k+1}-b_k}{\frac12(h_{k+1}+h_k)}.$$

The layer integration of the pressure gradient force is discussed in
detail by [@BURCHARDea97].

A conservative formulation can be derived for the recalculation of the
physical vertical velocity $w$ which is convenient in the discrete space
if $w$ is evaluated at the layer centres (see [@DELEERSNIJDERea92]):

$$\label{conservative_w}
w_k=\frac{1}{h_k}\left(
\partial_t(h_kz_{k-1/2})+\partial_x(p_kz_{k-1/2})+\partial_y(q_kz_{k-1/2})
+\bar w_kz_k-\bar w_{k-1}z_{k-1}\right).$$

It should be mentioned that $w$ only needs to be evaluated for
post-processing reasons.

For the layer-integrated tracer concentrations, we obtain the following
expression:

$$\label{C_Layer_Int}
\begin{array}{l}
\partial_t (h_k c^i_k) + \partial_x (p_kc^i_k)+\partial_y (q_k c^i_k)
+(\bar w_k+w^s_k) \tilde c^i_k 
- (\bar w_{k-1}+w^s_{k-1}) \tilde c^i_{k-1}\\
\\ \qquad
-(\nu'_t\partial_z c^i)_{k}
+(\nu'_t\partial_z c^i)_{k-1}
-\partial_x\left(A_k^Th_k\partial_xc^i_k\right)
-\partial_y\left(A_k^Th_k\partial_yc^i_k\right)
=Q^i_k.
\end{array}$$

It should be noted that the “horizontal” diffusion does no longer occur
along geopotential surfaces but along horizontal coordinate lines. The
properly transformed formulation would include some cross-diagonal terms
which may lead to numerical instabilities due to violation of
monotonicity. For an in-depth discussion of this problem, see
[@BECKERSea98] and [@BECKERSea00].

Horizontal curvilinear coordinates {#SectionCurviCoords}
----------------------------------

In this section, the layer-integrated equations from section
\[SectionTransformations\] are transformed to horizontal orthogonal
curvilinear coordinates. Similarly to general coordinates in the
vertical, these allow for much more flexibility when optimising
horizontal grids to coast-lines and bathymetry. Furthermore, this type
of coordinates system includes spherical coordinates as a special case.
The derivation of the transformed equations is carried out here
according to [@HAIDVOGELea99], see also [@ARAKAWAea77].

A rectangular domain with non-dimensional side lengths and with local
Cartesian coordinates ${\cal X}$ and ${\cal Y}$ is mapped to a physical
domain with four corners in such a way that the local coordinates of the
physical space, $(\xi_x,\xi_y)$ are orthogonal to each others
everywhere:

$$\label{calXYdef} 
{\cal X} \rightarrow \xi_x,\quad {\cal Y} \rightarrow \xi_y.$$

The infinitesimal increments in the physical space, $d\,\xi_x$ and
$d\,\xi_y$ are related to the infinitesimal increments in the
transformed space, $d\,{\cal X}$ and $d\,{\cal Y}$ by so-called metric
coefficients $m(x,y)$ and $n(x,y)$:

$$\label{curvidef} 
d\,\xi_x = \left(\frac{1}{m} \right) d\,{\cal X}, \quad
d\,\xi_y = \left(\frac{1}{n} \right) d\,{\cal Y}.$$

These metric coefficients have the physical unit of \[m$^{-1}$\]. With
$m=n=$const, Cartesian coordinates are retained, and with

$$\label{rE} 
m=\frac{1}{r_E\cos\phi},\quad n=\frac{1}{r_E},$$

spherical coordinates with ${\cal X}=\lambda$\[lambda\] and
${\cal Y}=\phi$ are retained (with the Earth’s radius $r_E$, longitude
$\lambda$ and latitude $\phi$).

With these notations, the layer-integrated equations
(\[ContiLayerInt\]), (\[uEqvi\]), and (\[vEqvi\]) given in section
\[SectionTransformations\] can be formulated as follows:\
[**Continuity equation:**]{}

$$\label{ContiLayerIntCurvi}
\partial_t \left(\frac{h_k}{mn}\right) 
+ \partial_{\cal X} \left(\frac{p_k}{n} \right)
+ \partial_{\cal Y} \left(\frac{q_k}{m} \right)
+ \frac{\bar w_k - \bar w_{k-1}}{mn}=0.$$

[**Momentum in $\xi_x$ direction:**]{} $$\label{uEqviCurvi}
\begin{array}{l}
\displaystyle 
\partial_t \left(\frac{p_k}{mn}\right)  
+\frac{\bar w_k \tilde u_k -\bar w_{k-1} \tilde u_{k-1}}{mn}  
-\frac{\tau^{\cal X}_k - \tau^{\cal X}_{k-1}}{mn}  
\\ \\ \quad
\displaystyle 
+\alpha\Bigg\{\partial_{\cal X}\left(\frac{u_kp_k}{n}\right)
+\partial_{\cal Y}\left(\frac{v_kp_k}{m}\right)
- q_k \left(\frac{f}{mn}+v_k\partial_{\cal X}\left(\frac{1}{n}\right)
-u_k\partial_{\cal Y}\left(\frac{1}{m}\right) \right)  
\\ \\ \quad
\displaystyle
-\partial_{\cal X}\left(\frac{2A_k^Mh_k}{n} m\partial_{\cal X}u_k\right)
-\partial_{\cal Y}\left(\frac{A_k^Mh_k}{m}
(n\partial_{\cal Y}u_k+m\partial_{\cal X}v_k)\right)
\\ \\ \quad
\displaystyle
-\frac{h_k}{n}\left(
\frac12h_N(\partial^*_{\cal X}b)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_{\cal X}b)_j
\right)\Bigg\}
=
-g\frac{h_k}{n}\partial_{\cal X}\zeta.
\end{array}$$

[**Momentum in $\xi_y$ direction:**]{} $$\label{vEqviCurvi}
\begin{array}{l}
\displaystyle 
\partial_t \left(\frac{q_k}{mn}\right)  
+\frac{\bar w_k \tilde v_k -\bar w_{k-1} \tilde v_{k-1}}{mn}  
-\frac{\tau^{\cal Y}_k - \tau^{\cal Y}_{k-1}}{mn}  
\\ \\ \quad
\displaystyle 
+\alpha\Bigg\{\partial_{\cal X}\left(\frac{u_kq_k}{n}\right)
+\partial_{\cal Y}\left(\frac{v_kq_k}{m}\right)
+ p_k \left(\frac{f}{mn}+v_k\partial_{\cal X}\left(\frac{1}{n}\right)
-u_k\partial_{\cal Y}\left(\frac{1}{m}\right) \right)  
\\ \\ \quad
\displaystyle
-\partial_{\cal Y}\left(\frac{2A_k^Mh_k}{m} n\partial_{\cal Y}v_k\right)
-\partial_{\cal X}\left(\frac{A_k^Mh_k}{n}
(n\partial_{\cal Y}u_k+m\partial_{\cal X}v_k)\right)
\\ \\ \quad
\displaystyle
-\frac{h_k}{m}\left(
\frac12h_N(\partial^*_{\cal Y}b)_N
+\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_{\cal Y}b)_j
\right)\Bigg\}
=
-g\frac{h_k}{m}\partial_{\cal Y}\zeta.
\end{array}$$

In (\[uEqviCurvi\]) and (\[vEqviCurvi\]), the velocity and momentum
components $u_k$ and $p_k$ are now pointing into the $\xi_x$-direction
and $v_k$ and $q_k$ are pointing into the $\xi_y$-direction. The
stresses $\tau^{\cal X}_k$ and $\tau^{\cal Y}_k$ are related to these
directions as well. In order to account for this rotation of the
velocity and momentum vectors, the rotational terms due to the Coriolis
rotation are extended by terms related to the gradients of the metric
coefficients. This rotation is here not considered for the horizontal
diffusion terms in order not to unnecessarily complicate the equations.
Instead we use the simplified formulation by [@KANTHAea00b], who argue
that it does not make sense to use complex formulations for minor
processes with highly empirical parameterisations.

Finally, the tracer equation is of the following form after the
transformation to curvilinear coordinates: $$\label{C_Layer_IntCurvi}
\begin{array}{l}
\displaystyle
\partial_t \left(\frac{h_k c^i_k}{mn}\right) 
+ \partial_{\cal X} \left(\frac{p_kc^i_k}{n}\right)
+\partial_{\cal Y} \left(\frac{q_k c^i_k}{m}\right)
+\frac{\bar w_k \tilde c^i_k - \bar w_{k-1} \tilde c^i_{k-1}}{mn}\\
\\ \qquad
\displaystyle
-\frac{(\nu'_t\partial_z c^i)_{k}
-(\nu'_t\partial_z c^i)_{k-1}}{mn}
\\ \\ \qquad
\displaystyle
-\partial_{\cal X}\left(\frac{A_k^Th_k}{n} m\partial_{\cal X}c^i_k\right)
-\partial_{\cal Y}\left(\frac{A_k^Th_k}{m} n\partial_{\cal Y}c^i_k\right)
=\frac{Q^i_k}{mn}.
\end{array}$$

Discretisation {#Section_discretisation}
==============

Mode splitting {#Section_mode_splitting}
--------------

The external system consisting of the surface elevation equation
(\[Elevation\]) and the transport equations (\[UMOM\]) and (\[VMOM\])
underlies a strict time step constraint if the discretisation is carried
out explicitly:

$$\label{CourantSurf}
\Delta t < \left[\frac12 \left(\frac{1}{\Delta x}+\frac{1}{\Delta y}\right)
\sqrt{2gD}\right]^{-1}.$$

In contrast to that, the time step of the internal system is only
depending on the Courant number for advection,

$$\label{CourantVel}
\Delta t < \min\left\{\frac{\Delta x}{u_{\max}},\frac{\Delta y}{v_{\max}}
\right\},$$

which in the case of sub-critical flow is a much weaker constraint. In
order not to punish the whole model with a small time step resulting
from the external system, two different approaches of mode splitting
have been developed in the past.

The first approach, in which the external mode is calculated implicitly,
has been proposed by [@MADALAea77]. This method is numerically stable
(if advection is absent) for unconditionally long time steps. The
temporal approximation is of second order if semi-implicit treatment is
chosen. In such models, the external and internal mode are generally
calculated with the same time steps (see e.g. [@BACKHAUS85]). The
introduction of interactions terms like (\[Slowfirst\]) - (\[Slowlast\])
is thus not necessary in such models.

Another approach is to use different time steps for the internal (macro
times steps $\Delta t$) and the external mode \[dtm\] (micro time steps
$\Delta t_m$). One of the first free surface models which has adopted
this method is the Princeton Ocean Model (POM), see [@BLUMBERGea87].
This method has the disadvantage that interaction terms are needed for
the external mode and that the consistency between internal and external
mode is difficult to obtain. The advantage of this method is that the
free surface elevation is temporally well resolved which is a major
requirement for models including flooding and drying. That is the reason
why this method is adopted here.

The micro time step $\Delta t_m$ has to be an integer fraction $M$ of
the macro time step $\Delta t$. $\Delta t_m$ is limited by the speed of
the surface waves (\[CourantSurf\]), $\Delta t$ is limited by the
current speed (\[CourantVel\]). The time stepping principle is shown in
figure \[figtimegrid\]. The vertically integrated transports are
averaged over each macro time step:

$$\label{Mdef} 
\bar U_{i,j}^{n+1/2} = \frac{1}{M}\sum_{l=n+0.5/M}^{n+(M-0.5)/M} U^l_{i,j}$$

and
$$\bar V_{i,j}^{n+1/2} = \frac{1}{M}\sum_{l=n+0.5/M}^{n+(M-0.5)/M} V^l_{i,j}$$
such that $$\begin{array}{l}
\displaystyle
\frac{\zeta_{i,j}^{n+1}-\zeta_{i,j}^{n}}{\Delta t}=
-\frac{\bar U_{i,j}^{n+1/2}-\bar U_{i-1,j}^{n+1/2}}{\Delta x}
-\frac{\bar V_{i,j}^{n+1/2}-\bar V_{i,j-1}^{n+1/2}}{\Delta y}.
\end{array}$$

![Sketch explaining the organisation of the time stepping.
[]{data-label="figtimegrid"}](./figures/gridtime.ps){width="12cm"}

Spatial discretisation {#Section_spatial_discretisation}
----------------------

For the spatial discretisation, a staggered C-grid is used, see
[@ARAKAWAea77]. The grid consists of prism-shaped finite volumes with
the edges aligned with coordinates. The reference grid for the tracer
points (from now on denoted by T-points) is shown in figures
\[fighorgrid\] and \[figvertgrid\]. The velocity points are located such
that the corresponding velocity components are centralised on the
surfaces of the T-point reference box, the $u$-velocity points (from now
on U-points) at the western and eastern surfaces, the $v$-velocity
points (from now on V-points) at the southern and northern surfaces and
the $w$-velocity points (from now on W-points) at the lower and upper
surfaces. The indexing is carried out with $i$-indices\[indexi\] in
eastern (${\cal X}$-) direction, with $j$-indices\[indexj\] in northern
(${\cal Y}$-) direction and with $k$-indices in upward ($z$-) direction,
such that each grid point is identified by a triple $(i,j,k)$. A T-point
and the corresponding eastern U-point, the northern V-point and the
above W-point have always the same index, see figures \[fighorgrid\] and
\[figvertgrid\]. The different grid points cover the following index
ranges:

$$\label{imaxjmax}
\begin{array}{llll}
\mbox{T-points:} & 1 \leq i \leq i_{\max}, & 
                   1 \leq j \leq j_{\max}, &
		   1 \leq k \leq k_{\max} \\
\mbox{U-points:} & 0 \leq i \leq i_{\max}, & 
                   1 \leq j \leq j_{\max}, &
		   1 \leq k \leq k_{\max} \\
\mbox{V-points:} & 1 \leq i \leq i_{\max}, & 
                   0 \leq j \leq j_{\max}, &
		   1 \leq k \leq k_{\max} \\
\mbox{W-points:} & 1 \leq i \leq i_{\max}, & 
                   1 \leq j \leq j_{\max}, &
		   0 \leq k \leq k_{\max} 
\end{array}$$

On the T-points, all tracers such as\[TSfirst\] temperature $T$,
salinity $S$, the general tracers $c^i$ and the density are located. All
turbulent quantities such as eddy viscosity $\nu_t$ and eddy diffusivity
$\nu_t'$ are located on the W-points.

\[lc\]\[\]\[1.4\][$(x_{i-1,j},y_{i-1,j})$]{}
\[lc\]\[\]\[1.4\][$(x_{i,j},y_{i,j})$]{}
\[lc\]\[\]\[1.4\][$(x_{i,j-1},y_{i,j-1})$]{}
\[lc\]\[\]\[1.4\][$(x_{i-1,j-1},y_{i-1,j-1})$]{} ![Grid layout and
indexing of corner points for curvilinear grids.
[]{data-label="fighorgridCurvi"}](./figures/CurviDXDY.ps "fig:"){width="8cm"}

For curvilinear grids, several arrays for spatial increments $\Delta x$
and $\Delta y$ have to be defined:

$$\label{dxdycuv+}
\begin{array}{rcl}
\Delta x^c_{i,j}&=&\left|\left|
\frac12(X_{i,j-1}+X_{i,j}-X_{i-1,j-1}-X_{i-1,j}) 
\right|\right|\\ \\
\Delta x^u_{i,j}&=&\left|\left|
\frac14(X_{i+1,j-1}+X_{i+1,j}-X_{i-1,j-1}-X_{i-1,j})\right|\right| \\ \\
\Delta x^v_{i,j}&=& \left|\left|X_{i,j}-X_{i-1,j}\right|\right| \\ \\
\Delta x^+_{i,j}&=&\left|\left|
\frac12(X_{i+1,j}-X_{i-1,j})\right|\right| \\ \\
\Delta y^c_{i,j}&=&\left|\left|
\frac12(X_{i-1,j}+X_{i,j}-X_{i-1,j-1}-X_{i,j-1})\right|\right| \\ \\
\Delta y^u_{i,j}&=&\left|\left|X_{i,j}-X_{i,j-1}\right|\right|
 \\ \\
\Delta y^v_{i,j}&=&\left|\left|
\frac14(X_{i-1,j+1}+X_{i,j+1}-X_{i-1,j-1}-X_{i,j-1})\right|\right| \\ \\
\Delta y^+_{i,j}&=&\left|\left|
\frac12(X_{i,j+1}-X_{i,j-1})\right|\right| \\ \\
\end{array}$$

where $\left|\left|X_{i,j}-X_{i-1,j}\right|\right|
=\left((x_{i,j}-x_{i-1,j})^2+(y_{i,j}-y_{i-1,j})^2\right) ^{1/2}$. The
superscripts $c,u,v,+$ in (\[dxdycuv+\]) indicate whether a $\Delta x$
or $\Delta y$ is centrered at a T-, U-, V-, or X-point, respectively.
For the locations of the corner points $X_{i,j}=(x_{i,j},y_{i,j})$, see
figure \[fighorgridCurvi\].

Lateral boundary conditions
---------------------------

Usually, a land mask is defined on the horizontal numerical grid. This
mask is denoted by $a^z$\[az\] for T-points, $a^u$\[au\] for U-points
and $a^v$\[av\] for V-points with $a^z$, $a^u$, and $a^v$ being integer
fields. A T-point is either a land point ($a^z=0$) or a water point
($a^z>0$). All U- and V-points surrounding a land point are defined as
closed boundary and masked out: $a^u=0$ and $a^v=0$. The velocities on
such closed boundaries are always set to 0.

Open boundaries are defined by $a^z>1$ for T-points. Forced boundary
points are marked by $a^z=2$ and passive boundary points by $a^z=3$. All
other T-points are characterised by $a^z=1$. For velocity points, three
different types are defined at the open boundaries. U-points are
classified by $a^u=3$ if both the T-points east and west are open
boundary points and by $a^u=2$ if one adjacent T-point is an open
boundary point and the other an open water point with $a^z=1$. The same
is carried out for V-points: They are classified by $a^v=3$ if both the
T-points south and north are open boundary points and by $a^v=2$ if one
adjacent T-point is an open boundary point and the other an open water
point with $a^z=1$. U-points which are adjacent to T-points with $a^z=2$
and which are not denoted by $a^u=2$ or $a^u=3$ are the external
U-points and are denoted by $a^u=4$. The same holds for V-points: Those
which are adjacent to T-points with $a^z=2$ and which are not denoted by
$a^v=2$ or $a^v=3$ are the external V-points and are denoted by $a^v=4$.

For a simple example of grid point classification, see figure \[mask\].

When the barotropic boundary forcing is carried out by means of
prescribed surface elevations only, then the surface elevation $\zeta$
is prescribed in all T-points with $a^z=2$. For passive boundary
conditions ($a^z=3$), where the curvature of the surface elevation is
zero normal to the boundary, the surface slope is simply extrapolated to
the boundary points. For a boundary point $(i,j)$ at the western
boundary this results e.g. in the following calculation for the boundary
point:

$$\zeta_{i,j}=\zeta_{i+1,j}+(\zeta_{i+1,j}-\zeta_{i+2,j})=
2\zeta_{i+1,j}-\zeta_{i+2,j}.$$

Bed friction {#SectionBedFric}
------------

As already mentioned earlier in section \[SectionDynBounds\], caution is
needed when discretising the bottom boundary conditions for momentum,
(\[DynBoundbx\]). They are an example for a physical condition which has
to be modified for the numerical discretisation, since the discrete
velocity point nearest to the bottom is half a grid box away from the
point where the boundary condition is defined. Furthermore, due to the
logarithmic law, high velocity gradients are typical near the bed.
Simply setting the discrete bottom velocity to zero, would therefore
lead to large discretisation errors. Instead, a flux condition using
bottom stresses is derived from the law of the wall.

For the determination of the normalised bottom stresses

$$\label{tauxb}
\frac{\tau^x_b}{\rho_0}=u_*^{bx}u_*^b,$$

$$\label{tauyb}
\frac{\tau^y_b}{\rho_0}=u_*^{by}u_*^b$$

with the friction velocities $u_*^b=\sqrt{\tau_b/\rho_0}$\[pagetaub\]
with $\tau_b=\sqrt{(\tau^x_{b})^2+(\tau^y_{b})^2}$, assumptions about
the structure of velocity inside the discrete bottom layer have to be
made. We use here the logarithmic profile

$$\label{log_prof}
\frac{u(z')}{u_*}
=\frac{1}{\kappa}\mbox{ln}\left(\frac{z'+z_0^b}{z_0^b}\right),$$

with the bottom roughness length $z_0^b$, the von Kármán constant
$\kappa=0.4$ and the distance from the bed, $z'$. Therefore, estimates
for the velocities in the centre of the bottom layer can be achieved by:

$$\label{ulogdis}
u_b = \frac{u_*^{bx}}{\kappa}\mbox{ln} \left(\frac{0.5h_1+z_0^b}{z_0^b}\right),$$

$$\label{vlogdis}
v_b = \frac{u_*^{by}}{\kappa}\mbox{ln} \left(\frac{0.5h_1+z_0^b}{z_0^b}\right).$$

For $h_1\rightarrow 0$, the original Dirichlet-type no-slip boundary
conditions (\[DynBoundbx\]) are retained. Another possibility would be
to specify the bottom velocities $u_b$ and $v_b$ such that they are
equal to the layer-averaged log-law velocities (see [@BAUMERTea92]). The
calculation of this is however slightly more time consuming and does not
lead to a higher accuracy.

Drying and flooding {#Section_dry}
-------------------

The main requirement for drying and flooding is that the vertically
integrated fluxes $U$ and $V$ are controlled such that at no point a
negative water depth occurs. It is clear that parts of the physics which
play an important role in very shallow water of a few centimetres depth
like non-hydrostatic effects are not included in the equations. However,
the model is designed in a way that the control of $U$ and $V$ in very
shallow water is mainly motivated by the physics included in the
equations rather than by defining complex drying and flooding
algorithms. It is assumed that the major process in this situation is a
balance between pressure gradient and bottom friction. Therefore, in the
case of very shallow water, all other terms are multiplied with the
non-dimensional factor $\alpha$ which approaches zero when a minimum
water depth is reached.

By using formulation (\[bottom\_vert\]) for calculating the bottom drag
coefficient $R$, it is guaranteed that $R$ is exponentially growing if
the water depth approaches very small values. This slows the flow down
when the water depth in a velocity point is sinking and also allows for
flooding without further manipulation.

In this context, one important question is how to calculated the depth
in the velocity points, $H^u$ and $H^v$, since this determines how
shallow the water in the velocity points may become on sloping beaches.
In ocean models, usually, the depth in the velocity points is calculated
as the mean of depths in adjacent elevation points (T-points):

$$\label{Hmean}
H^u_{i,j}=\frac12\left(H_{i,j}+H_{i+1,j}\right),
\qquad
H^v_{i,j}=\frac12\left(H_{i,j}+H_{i,j+1}\right).$$

Other models which deal with drying and flooding such as the models of
[@DUWE88] and [@CASULLIea94] use the minimum of the adjacent depths in
the T-points:

$$\label{Hmin}
H^u_{i,j}= \min\{H_{i,j},H_{i+1,j}\},
\qquad
H^v_{i,j}= \min\{H_{i,j},H_{i,j+1}\}.$$

This guarantees that all depths in the velocity points around a T-point
are not deeper than the depth in the T-point. Thus, when the T-point
depth is approaching the minimum depth, then all depths in the velocity
points are also small and the friction coefficient correspondingly
large.

Each of the methods has however drawbacks: When the mean is taken as in
(\[Hmean\]), the risk of negative water depths is relatively big, and
thus higher values of $D_{\min}$ have to be chosen. When the minimum is
taken, large mud-flats might need unrealistically long times for drying
since all the water volume has to flow through relatively shallow
velocity boxes. Also, velocities in these shallow boxes tend to be
relatively high in order to provide sufficient transports. This might
lead to numerical instabilities.

Therefore, GETM has both options, (\[Hmean\]) and (\[Hmin\]) and the
addition of various other options such as depth depending weighting of
the averaging can easily be added. These options are controlled by the
GETM variable [vel\_depth\_method]{}, see section \[sec-uv-depth\]
(subroutine [uv\_depths]{}) documented on page .

If a pressure point is dry (i.e. its bathymetry value is higher than a
neighbouring sea surface elevation), the pressure gradient would be
unnaturally high with the consequence of unwanted flow acceleration.
Therefore this pressure gradient will be manipulated such that (only for
the pressure gradient calculation) a virtual sea surface elevation
$\tilde \zeta$ is assumed (see figure \[figpressgrad\]). In the
situation shown in figure \[figpressgrad\], the left pressure point is
dry, and the sea surface elevation there is for numerical reasons even
slightly below the critical value $-H_{i,j}+H_{\min}$. In order not to
let more water flow out of the left cell, the pressure gradient between
the two boxes shown is calculated with a manipulated sea surface
elevation on the right, $\tilde \zeta_{i+1,j}$.

See also [@BURCHARDea03a] for a description of drying and flooding
numerics in GETM.

Introduction to the calculation domain
======================================

This module handles all tasks related to the definition of the
computational domain - except reading in variables from file. The
required information depends on the $grid\_type$ and also on the
complexity of the model simulation to be done. The mandatory varible
$grid\_type$ read from the file containing the bathymetry and coordinate
information (presently only NetCDF is supported) is guiding subsequent
tasks. $grid\_type$ can take the following values:

-   equi-distant plane grid - $dx$, $dy$ are constant - but not
    necessarily equal

-   equi-distant spherical grid - $dlon$, $dlat$ are constant - and
    again not necessarily equal

-   curvilinear grid in the plane - $dx$, $dy$ are both functions of
    (i,j). The grid must be orthogonal

For all values of $grid\_type$ the bathymetry given on the T-points (see
the GETM manual for definition) must be given.

Based on the value of $grid\_type$ the following additional variables
are required:

-   proper monotone coordinate informtion in the xy-plane with
    equidistant spacing. The name of the coordinate variables are
    $xcord$ and $ycord$.

-   proper monotone coordinate informtion on the sphere with equidistant
    spacing in longitude and latitude. The names of the coordinate
    variables are $xcord$ and $ycord$.

-   position in the plane of the grid-vertices. These are called
    X-points in GETM. The names of these two variables are $xx$ and
    $yx$.

In addition to the above required grid information the following
information is necessary for specific model configurations:

-   $latu$ and $latv$ If $f\_plane$ is false information about the
    latitude of U- and V-points are required for calculating the
    Coriolis term correctly. For $grid\_type=1$ $latu$ and $latv$ are
    calculated based on an additional field $latc$ i.e. the latitude of
    the T-points. For $grid\_type=3$ $latx$ i.e. the latitude of the
    X-points will have to be provided in order to calculate $latu$ and
    $latv$.

-   $lonc$, $latc$ and $convc$ The longitude, latitude positions of the
    T-points are required when using forcing from a NWP-model. $lonc$
    and $latc$ are used to do spatial interpolation from the meteo-grid
    to the GETM model and $convc$ is the rotation of the local grid from
    true north.

In addition to the information above a few files are optionally read in
$init\_domain()$. Information about open boundaries, modifications to
the bathymetry and the calculation masks are are done via simple ASCII
files.

Introduction to 2d module
=========================

In the 2D module of GETM the vertically integrated mode is calculated,
which is basically the vertically integrated momentum equations and the
sea surface elevation equation. For the momentum equations, interaction
terms with the baroclinic three-dimentional mode need to be considered.
Those terms are here called the slow terms.

Vertically integrated mode {#SectionVerticalIntegrated}
--------------------------

In order to provide the splitting of the model into an internal and an
external mode, the continuity equation and the momentum equations are
vertically integrated. The vertical integral of the continuity equation
together with the kinematic boundary conditions (\[KinBoundsSurf\]) and
(\[KinBoundsBott\]) gives the sea surface elevation equation:

$$\label{Elevation}
\partial_t \zeta = - \partial_x U- \partial_y V.$$

with

$$U=\int_{-H}^{\zeta} u\,dz,\qquad V=\int_{-H}^{\zeta} v\,dz.$$

Integrating the momentum equations (\[uEq\]) and (\[vEq\]) vertically
results in:

$$\label{Uint}
\begin{array}{l}
\displaystyle
\partial_tU+
\tau_b^x+\alpha\bigg(\int_{-H}^{\zeta}
\left(\partial_x u^2+\partial_y(uv)\right)\,dz\\
\\
\displaystyle
\quad-\tau^x_s
-\int_{-H}^{\zeta}\big(
\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle
\quad -fV
-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_x b\,dz'\,dz \bigg)
=
-gD\partial_x\zeta
\end{array}$$

and

$$\label{Vint}
\begin{array}{l}
\displaystyle
\partial_tV+
\tau_b^y+\alpha\bigg(\int_{-H}^{\zeta}
\left(\partial_x(uv)+\partial_yv^2\right))\,dz \\
\\
\displaystyle
\quad-\tau^y_s
-\int_{-H}^{\zeta}\big(
\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle
\quad +fU
-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_y b\,dz'\,dz \bigg)
=
-gD\partial_y\zeta.
\end{array}$$

Here, $\tau_b^x$ and $\tau_b^y$ are bottom stresses. Their calculation
is discussed in section \[sec-bottom-friction-3d\]. As a first
preparation for the mode splitting, these integrals of the momentum
equations can be formally rewritten as

$$\label{UMOM}
\begin{array}{l}
\displaystyle
\partial_tU+\frac{R}{D^2}U\sqrt{U^2+V^2}+S^x_{F}
+\alpha\bigg(\partial_x\left(\frac{U^2}{D}\right)
+\partial_y\left(\frac{UV}{D}\right)\\
\\
\displaystyle
\quad-\tau^x_s
-\partial_x\left(2A_h^MD\partial_x\left(\frac{U}{D}\right)\right)
-\partial_y\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right)
\\ \\
\displaystyle
\quad-fV
+S^x_{A}-S^x_D+S^x_B\bigg)
=
-gD\partial_x\zeta
\end{array}$$

and

$$\label{VMOM}
\begin{array}{l}
\displaystyle
\partial_tV+\frac{R}{D^2}V\sqrt{U^2+V^2}+S^y_{F}
+\alpha\bigg(\partial_x\frac{UV}{D}+\partial_y\frac{V^2}{D}\\
\\
\displaystyle
\quad-\tau^y_s
-\partial_x\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right)
-\partial_y\left(2A_h^MD\partial_y\left(\frac{V}{D}\right)\right)
\\ \\
\displaystyle
\quad+fU
+S^y_{A}-S^y_D+S^y_{B}\bigg)
=
-gD\partial_y\zeta
\end{array}$$

with the so-called slow terms for bottom friction

$$\label{Slowfirst}
S^x_{F}=\tau^x_b-\frac{R}{D^2}U\sqrt{U^2+V^2},$$

$$\label{Slowsecond}
S^y_{F}=\tau^y_b-\frac{R}{D^2}V\sqrt{U^2+V^2},$$

horizontal advection

$$\label{SxA}
S^x_{A}=\int_{-H}^{\zeta}\left(\partial_x u^2+\partial_y(uv)\right)\,dz-
\partial_x\left(\frac{U^2}{D}\right)-\partial_y\left(\frac{UV}{D}\right),$$

$$\label{SyA} 
S^y_{A}=\int_{-H}^{\zeta}\left(\partial_x (uv)+\partial_yv^2\right)\,dz-
\partial_x\left(\frac{UV}{D}\right)-\partial_y\left(frac{V^2}{D}\right),$$

horizontal diffusion

$$\label{SxD}
\begin{array}{l}
\displaystyle
S^x_D=
\int_{-H}^{\zeta}\big(
\partial_x\left(2A_h^M\partial_xu\right)-\partial_y\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle\qquad 
-\partial_x\left(2A_h^MD\partial_x\left(\frac{U}{D}\right)\right)
-\partial_y\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right), 
\end{array}$$

$$\label{SyD}
\begin{array}{l}
\displaystyle
S^y_D=
\int_{-H}^{\zeta}\big(
\partial_y\left(2A_h^M\partial_yv\right)-\partial_x\left(A_h^M
(\partial_yu+\partial_xv)\right)\big)\,dz
\\ \\
\displaystyle\qquad 
-\partial_y\left(2A_h^MD\partial_y\left(\frac{V}{D}\right)\right)
-\partial_x\left(A_h^MD \left(\partial_y\left(\frac{U}{D}\right)
+\partial_x\left(\frac{V}{D}\right)\right)\right), 
\end{array}$$

and internal pressure gradients

$$S^x_B=-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_x b\,dz'\,dz$$

and

$$\label{Slowlast}
S^y_B=-\int_{-H}^{\zeta}\int_z^{\zeta}\partial_y b\,dz'\,dz.$$

The drag coefficient $R$ for the external mode is calculated as (this
logarithmic dependence of the bottom drag from the water depth and the
bottom roughness parameter $z_b^0$ is discussed in detail by
[@BURCHARDea02]):

$$\label{bottom_vert}
R = \left(\frac{\kappa}
{\ln\left(\frac{\frac{D}{2}+z^b_0}{z^b_0}\right)}\right)^2.$$

It should be noted that for numerical reasons, an additional explicit
damping has been implemented into GETM. This method is based on
diffusion of horizontal transports and is described in section
\[sec-uv-diffusion\] on page .

Introduction to 3d module
=========================

Overview over 3D routines in GETM {#Section_Overview_3D}
---------------------------------

This module contains the physical core of GETM. All three-dimensional
equations are iterated here, which are currently the equations for

  ---------- -------------------------- --------------- ---------- --------------------- ------
  quantity   description                unit            variable   routine name          page
  $p_k$      layer-int. $u$-transport   m$^2$s$^{-1}$   [uu]{}     [uu\_momentum]{}      
  $q_k$      layer-int. $v$-transport   m$^2$s$^{-1}$   [vv]{}     [vv\_momentum]{}      
  $\theta$   potential temperature      $^{\circ}$C     [T]{}      [do\_temperature]{}   
  $S$        salinity                   psu             [S]{}      [do\_salinity]{}      
  $C$        suspended matter           kgm$^{-3}$      [spm]{}    [do\_spm]{}           
  ---------- -------------------------- --------------- ---------- --------------------- ------

The vertical grid for GETM, i.e. the layer thicknesses in all U-, V- and
T-points, are defined in the routine [coordinates]{}, see section
\[sec-coordinates\] on page \[sec-coordinates\].

The grid-related vertical velocity $\bar w_k$ is calculated directly
from the layer-integrated continuity equation (\[ContiLayerInt\]) which
here done in the routine [ww\_momentum]{} described on page .

The physics of the horizontal momentum equations is given in section
\[Section\_3d\_momentum\], and their transformation to general vertical
coordinates in section \[SectionLayerIntegrated\]. Their numerical
treatment will be discussed in the routines for the individual terms,
see below. The forcing terms of the horizontal momentum equations are
calculated in various routines, such as [uv\_advect\_3d]{} for the
three-dimensional advection (which in turn calls [advection\_3d]{} in
case that higher order positive definite advection schemes are chosen
for the momentum equation), [uv\_diffusion\_3d.F90]{} for the horizontal
diffusion, [bottom\_friction\_3d]{} for the bottom friction applied to
the lowest layer, and [internal\_pressure]{} for the calculation of the
internal pressure gradients.

The major tracer equations in any ocean model are those for potential
temperature and salinity. They are calculated in the routines
[do\_temperature]{} and [do\_salinity]{}. A further hard-coded tracer
equation is the suspended matter equation, see [do\_spm]{}.

In the near future (the present text is typed in February 2006), a
general interface to the biogeochemical module of GOTM (also not yet
released) will be available. This allow to add tracer equations of
arbitrary complexity to GETM, ranging from completely passive tracer
equations to complex ecosystem models such as ERSEM ([@BARETTAea95]).
The interfacing between this so-called GOTM-BIO to GETM is made in a
similar manner than the interfacing between GETM and the GOTM turbulence
module described in [gotm]{} on page . The basic structure of GOTM-BIO
has been recently presented by [@BURCHARDea06]. Some more details about
the tracer equations currently included in GETM is given in section
\[Section\_tracer\].

The entire turbulence model, which basically provides eddy viscosity
$\nu_t$ and eddy diffusivity $\nu'_t$ is provided from the General Ocean
Turbulence Model (GOTM, see [@UMLAUFea05] for the source code
documentation and [http://www.gotm.net]{} download of source code,
docomentation and test scenarios). The turbulence module of GOTM (which
is a complete one-dimensional water column model) is coupled to GETM via
the interfacing routine [gotm]{} described in section [gotm]{} on page .
Major input to the turbulence model are the shear squared
$M^2=\left(\partial_zu\right)^2+\left(\partial_zu\right)^2$ and the
buoyancy frequency squared $N^2=\partial_z b$ with the buoyancy $b$ from
(\[bdef\]). Those are calculated and interpolated to the T-points where
the turbulence model columns are located in the routine [ss\_nn]{}
described on page .

The surface and bottom stresses which need to be passed to the
turbulence module as well, are interpolated to T-points in the routine
[stresses\_3d]{}, see page .

The module [rivers]{} (see section \[sec-rivers\] on page ) organises
the riverine input of fresh water from any number of rivers.

Three-dimensional boundary conditions for temperature and salinity are
provided by means of the module [bdy-3d]{}, see section \[bdy-3d\]
described on page .

The remaining routines in the module [3d]{} deal with the coupling of
the external and the internal mode. The basic idea of the mode splitting
has already been discussed in section \[Section\_mode\_splitting\]. The
consistency of the two modes is given through the so-called slow terms,
which are mode interaction terms resulting from subtracting vertically
integrated equations with parameterised advection, diffusion, bottom
friction and internal pressure gradient from vertically integrated
equations with explicit vertical resolution of these processes. These
slow terms which are updated every macro time step only (that is why we
call them slow terms) need to be considered for the external mode
included in module [2d]{}. Those slow terms are calculated here in the
[3d]{} module at the end of [integrate\_3d]{} and in the routine
[slow\_bottom\_friction]{}, and they are added together in
[slow\_terms]{}, see the descriptions in sections \[sec-integrate-3d\],
\[sec-slow-bottom-friction\] and \[sec-slow-terms\] on pages , , and ,
respectively.

One other important measure of coupling the two modes is to add to all
calculated $u$- and $v$-velocity profiles the difference between their
vertical integral and the time-average of the vertically integrated
transport from the previous set of micro time steps. This shifting is
done in the routines [uu\_momentum\_3d]{} and [vv\_momentum\_3d]{} and
the time-average of the vertically integrated transport is updated in
the [2d]{} module in the routine [m2d]{} and divided by the number of
micro time steps per macro time step in [start\_macro]{}. Further basic
calculations performed in [start\_macro]{} (see description in section
\[sec-start-macro\] on page ) are the updates of the [*old*]{} and
[*new*]{} sea surface elevations with respect to the actual macro time
step. The routine [stop\_macro]{} (see description in section
\[sec-stop-macro\] on page ) which called at the end of each macro time
step simply resets the variables for the time-averaged transports to
zero.

Tracer equations {#Section_tracer}
----------------

The general conservation equation for tracers $c^i$ with
$1\leq i \leq N_c$ (with $N_c$ being the number of tracers), which can
e.g. be temperature, salinity, nutrients, phytoplankton, zoo-plankton,
suspended matter, chemical concentrations etc. is given as:

$$\label{densz}
\begin{array}{l}
\partial_t c^i +\partial_x (uc^i) +\partial_y(vc^i) 
+\partial_z ((w+\alpha w_s^i)c^i)
-\partial_z(\nu_t' \partial_z c^i)
\\ \\ \displaystyle \qquad
-\partial_x(A_h^T \partial_x c^i)
-\partial_y(A_h^T \partial_y c^i)
=Q^i.
\end{array}$$

Here, $\nu'_t$ denotes the vertical eddy diffusivity and $A_h^T$ the
horizontal diffusivity. Vertical migration of concentration with
migration velocity $w_s^i$ (positive for upward motion) is considered as
well. This could be i.e. settling of suspended matter or active
migration of phytoplankton. In order to avoid stability problems with
vertical advection when intertidal flats are drying, the settling of SPM
is linearly reduced towards zero when the water depth is between the
critical and the minimum water depth. This is done by means of
multiplication of the settling velocity with $\alpha$, (see the
definition in equation (\[alpha\])). $Q^i$ denotes all internal sources
and sinks of the tracer $c^i$. This might e.g. be for the temperature
equation the heating of water due to absorption of solar radiation in
the water column.

Surface of bottom boundary conditions for tracers are usually given by
prescribed fluxes:

$$\label{SurfFlux}
-\alpha w_s^i c^i+\nu_t' \partial_z c^i = F^i_s  \qquad \mbox{for } z=\zeta$$

and

$$\label{BottFlux}
-\alpha w_s^i c^i+\nu_t' \partial_z c^i = -F^i_b  \qquad \mbox{for } z=-H,$$

with surface and bottom fluxes $F^n_s$ and $F^n_b$ directed into the
domain, respectively.

At open lateral boundaries, the tracers $c^n$ are prescribed for the
horizontal velocity normal to the open boundary flowing into the domain.
In case of outflow, a zero-gradient condition is used.

All tracer equations except those for temperature, salinity and
suspended matter will be treated in the future by means of GOTM-BIO.

The two most important tracer equations which are hard-coded in GETM are
the transport equations for potential temperature $\theta$ in
$^{\circ}$C and salinity $S$ in psu (practical salinity units):

$$\label{Teq}
\begin{array}{l}
\partial_t \theta +\partial_x (u\theta) +\partial_y(v\theta) +\partial_z (w\theta)
-\partial_z(\nu_t' \partial_z \theta)
\\ \\ \displaystyle \qquad
-\partial_x(A_h^\theta \partial_x \theta)
-\partial_y(A_h^\theta \partial_y \theta)
=\frac{\partial_z I}{c'_p \rho_0},
\end{array}$$

$$\label{Seq}
\begin{array}{l}
\partial_t S +\partial_x (uS) +\partial_y(vS) +\partial_z (wS)
-\partial_z(\nu_t' \partial_z S)
\\ \\ \displaystyle \qquad
-\partial_x(A_h^S \partial_x S)
-\partial_y(A_h^S \partial_y S)
=0.
\end{array}$$

On the right hand side of the temperature equation (\[Teq\]) is a source
term for absorption of solar radiation with the solar radiation at depth
$z$, $I$, and the specific heat capacity of water, $c'_p$. According to
[@PAULSONea77] the radiation $I$ in the upper water column may be
parameterised by

$$\label{Light}
I(z) = I_0 \left(ae^{-\eta_1z}+(1-a)e^{-\eta_2z}\right).$$

Here, $I_0$ is the albedo corrected radiation normal to the sea surface.
The weighting parameter $a$ and the attenuation lengths for the longer
and the shorter fraction of the short-wave radiation, $\eta_1$ and
$\eta_2$, respectively, depend on the turbidity of the water.
[@JERLOV68] defined 6 different classes of water from which
[@PAULSONea77] calculated weighting parameter $a$ and attenuation
coefficients $\eta_1$ and $\eta_2$.

At the surface, flux boundary conditions for $T$ and $S$ have to be
prescribed. For the potential temperature, it is of the following form:

$$\label{TempFlux}
\nu'_t \partial_z T= \frac{Q_s+Q_l+Q_b}{c'_p \rho_0},
\qquad \mbox{for } z=\zeta,$$

with the sensible heat flux, $Q_s$, the latent heat flux, $Q_l$ and the
long wave back radiation, $Q_b$. Here, the [@KONDO75] bulk formulae have
been used for calculating the momentum and temperature surface fluxes
due to air-sea interactions. In the presence of sea ice, these air-sea
fluxes have to be considerably changed, see e.g. [@KANTHAea00b]. Since
there is no sea-ice model coupled to GETM presently, the surface heat
flux is limited to positive values, when the sea surface temperature
$T_s$ reaches the freezing point

$$\label{freezingpoint}
T_f=-0.0575\,S_s+1.710523\cdot 10^{-3}\, S_s^{1.5}
-2.154996\cdot 10^{-4}\,S_s^2\approx -0.0575\,S_s$$

with the sea surface salinity $S_s$, see e.g. [@KANTHAea00]:

$$\label{primitive_ice_model} 
Q_{surf} = \left\{
\begin{array}{ll}
Q_s+Q_l+Q_b, & \mbox{ for } T_s > T_f, \\ \\
\max\{0,Q_s+Q_l+Q_b\}, & \mbox{ else.}
\end{array}
\right.$$

For the surface freshwater flux, which defines the salinity flux, the
difference between evaporation $Q_E$ (from bulk formulae) and
precipitation $Q_P$ (from observations or atmospheric models) is
calculated:

$$\label{SalFlux}
\nu'_t\partial_z S = \frac{S(Q_E-Q_P)}{\rho_0(0)},
\qquad \mbox{for } z=\zeta,$$

where $\rho_0(0)$ is the density of freshwater at sea surface
temperature. In the presence of sea-ice, the calculation of freshwater
flux is more complex, see e.g. [@LARGEea94]. However, for many short
term calculations, the freshwater flux can often be neglected compared
to the surface heat flux.

A complete revision of the surface flux calculation is currently under
development. It will be the idea to have the same surface flux
calculations for GOTM and GETM. In addition to the older bulk formulae
by [@KONDO75] we will also implement the more recent formlations by
[@FAIRALLea96].

Heat and salinity fluxes at the bottom are set to zero.

Equation of state {#Section_state_eq}
-----------------

The coupling between the potential temperature and salinity equations
and the momentum equations is due to an algebraic equation of state:

$$\label{UNESCO} 
\rho=\rho(\theta,S,p_0)$$

with $p_0=g\rho_0(\zeta-z)$ being the hydrostatic reference pressure. In
order to obtain potential density from the equation of state, $p_0$
needs to be set to zero, which is the default in GETM.

Currently the equation of state by [@FOFONOFFea83] is implemented into
GETM, but the more recent and more consistent equation of state by
[@JACKETTea05] which is already contained in GOTM will be added as an
option in the near future.

For the equation of state, also linearised version are implemented into
GETM, for details, see section \[sec-eqstate\] on page .

For convinient use in other subroutines the buoyancy $b$ as defined in
is calculated and stored in the GETM variable [buoy]{}.

NetCDF I/O modules
==================

The use of external files - both input and output - is done via generic
wrapper routines in GETM. For specific formats the I/O routines must be
coded. In this section the specific NetCDF related I/O routines are
given.
