Here, the general vertical coordinates layer distribution
in T-, U- and V-points is calculated. The general vertical coordinates as
discussed in section \ref{SectionGeneralCoordinates},
see equations (\ref{sigma}) - (\ref{MLDTransform}), are basically
an interpolation between equidistant and non-equaidistant $\sigma$
coordinates. During the first call, a three-dimensional field
{\tt gga} containing the relative interface positions is calculated,
which further down used together with the actual water depth in the
T-, U- and V-points for calculating the updated old and new layer
thicknesses.

