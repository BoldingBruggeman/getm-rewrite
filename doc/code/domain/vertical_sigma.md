Here, the sigma coordinates layer distribution in T-, U- and V-points is calculated.
The layer interfaces for each
layer index have a fixed relative position $\sigma_k$ in the water column,
which may be even equidistant or non-equidistant, see equations
(\ref{sigma}) and (\ref{formula_Antoine}).
The surface and bottom zooming factors
$d_u$ and $d_l$ are read in via the {\tt domain} namelist in {\tt getm.inp}
as {\tt ddu} and {\tt ddl}.
In the first call to the {\tt sigma\_coordinates}, the relative interface positions
{\tt dga} are calculated as a one-dimensional vector (in case of
non-equidistant $\sigma$ coordinates), and those are then multiplied with
the water depths in all T-, U- and V-points to get the layer thicknesses.

