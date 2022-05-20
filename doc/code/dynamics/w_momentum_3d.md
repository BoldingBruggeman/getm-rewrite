Here, the local continuity equation is calculated in order to obtain
the grid-related vertical velocity $\bar w_k$. An layer-integrated equation
for this quantity is given as equation (\ref{ContiLayerInt}) which
has been derived
from the differential formulation (\ref{Konti}).

Since the kinematic boundary condition must hold (and is used for the
derivation of (\ref{ContiLayerInt})), the grid-related vertical
velocity at the surface muzst be zero, i.e.\ $\bar w_{k_{\max}}=0$.
This is a good consistence check for the mode splitting, since this is
only fulfilled if the vertically integrated continuity equation
(which is the sea surface elevation equation calculated on the
micro time step) and this local continuity equation are compatible.

The physical vertical velocity is then recalculated from the grid-related
vertical velocity by means of (\ref{conservative_w}), ... which should
soon be coded in the routine {\tt tow} in the directory {\tt futils}.

