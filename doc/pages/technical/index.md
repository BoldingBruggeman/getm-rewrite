title: Technical implementation
Author: Karste Bolding, Jorn Bruggeman and Hans Burchard

The technical implementation of GETM.

Mainly for Hans - maybe need to add ssh key

 1. cd
 2. mkdir -p source/repos && cd ~/source/repos
 3. git clone --recurse-submodules git@github.com:BoldingBruggeman/getm-rewrite.git
 4. cd getm-rewrite/
 6. git checkout devel
 7. mkdir _build && cd _build && cmake ..
 8. make


### Loop boundaries and halo updates

The calculation domain is defined by the indices _(imin:imax,jmin:jmax,kmin:kmax)_. The sizes of the numerical arrays are given by _(imin-halo:imax+halo,jmin-halo:jmax+halo,kmin:kmax)_ where _halo_ is the width of the socalled _halo-zones_ . The halo-zones are used when parallizing the sompuation using sub-domain decomposition.

Basically the solution of the equations are done by a series of do-loops over the model domain for the different components and equations involved. A concise implementation of loop boundaries help increase the effieciency of the implementation - and - implies fewer exchanges between the subdomains.

The model variables falls in four different groups:

  1. calculation domain related - defined on T, U, V and X grids
  2. state variables on U, V and X grids
  3. aaaa 
  4. auxilliary variables needed to help calculate/update state-variables

The first group of variables are defined in the source file [[domain.F90]] and the derived types [[type_getm_grid]] and [[type_getm_domain]]. Some of these variables are only calculated during initialization of the calculation domain.

| Variable         |     Loop boundaries                        | Halo updates |
|:----------------:|:------------------------------------------:|:------------:|
| \(H\)            |           none - read from file            |       8      |
| \(H_u,H_v,H_x\) ?|           imin:imax,jmin:jmax              |       8      |
| \(H_u\)     ??   | imin-halo:imax+halo-1,jmin-halo:jmax+halo  |       3      |
| \(H_v\)     ??   | imin-halo:imax+halo,jmin-halo:jmax+halo-1  |       3      |
| \(H_x\)     ??   | imin-halo:imax+halo-1,jmin-halo:jmax+halo-1|       8      |
| \(H_u,H_v\)      |           imin:imax,jmin:jmax              |       8      |
| \(z\)            |           imin:imax,jmin:jmax              |       8      |
| \(D,h_n\)        |  imin-halo:imax+halo,jmin-halo:jmax+halo   |       0      |
| \(z_u,D_u,h_u\)  | imin-halo:imax+halo-1,jmin-halo:jmax+halo  |       3      |
| \(z_v,D_v,h_v\)  | imin-halo:imax+halo,jmin-halo:jmax+halo-1  |       3      |
| \(D_x,h_x\)      |           imin:imax,jmin:jmax              |       8      |
|                  |                                            |              |
|                  |                                            |              |


The second group of variables are defined in the source file [[momentum.F90]] and the derived types [[type_getm_momentum]] and ????[[type_getm_domain]].

| Variable         |     Loop boundaries                       | Halo updates  |
|:----------------:|:-----------------------------------------:|:-------------:|
| \(U,V\)          |           imin:imax,jmin:jmax             |       8       |
| \(p_k,q_k\)      |           imin:imax,jmin:jmax:k           |       8       |
| \(u_1,v_1\)      | imin-halo:imax+halo,jmin-halo:jmax+halo   |       0       |
| \(u_k,v_k\)      | imin-halo:imax+halo,jmin-halo:jmax+halo:k |       0       |
|                  |                                           |               |
|                  |                                           |               |


### Code review

Domain

[[uv_depths.F90]] (open)

[[cfl_check.F90]] (open)

[[metrics.F90]] (open - only cartesian for now)

[[depth_update.F90]] (open)

are drying variables only needed inside calculation domain (imin:imax,jmin:jmax)?

[[test_uv_depths]]

[[test_uv_depths]]

[[test_advection]]

2D

[[sealevel.F90]] (open)

[[pressure_surface.F90]] (open)

[[momentum_2d.F90]] (open)

3D

[[temperature.F90]] (open - nothing done yet - only definitions)

[[salinity.F90]] (open)

[[radiation.F90]] (open)

[[density.F90]] (open)

[[mixing.F90]] (open)

### Advection implementation

[[operators.F90]] (open)

[[advection.F90]] (open)

[[operators.F90]] (open)

[[advection.F90]] (open)

[[advection.F90.template]] (open)

[[advection_superbee.F90]] (open)

