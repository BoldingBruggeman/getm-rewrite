! Size of dimensions with and without HALO-zones:
#define HALO 2
#define _IRANGE_HALO_ imin-HALO:imax+HALO
#define _JRANGE_HALO_ jmin-HALO:jmax+HALO
#define _IRANGE_NO_HALO_ imin:imax
#define _JRANGE_NO_HALO_ jmin:jmax
#define _KRANGE_ 0:kmax

#define E2DFIELD  _IRANGE_HALO_,_JRANGE_HALO_
#define I2DFIELD  _IRANGE_HALO_,_JRANGE_HALO_
#define I3DFIELD  _IRANGE_HALO_,_JRANGE_HALO_,_KRANGE_

#ifdef _WRITE_HALOS_
#define _IRANGE_ _IRANGE_HALO_
#define _JRANGE_ _JRANGE_HALO_
#else
#define _IRANGE_ _IRANGE_NO_HALO_
#define _JRANGE_ _JRANGE_NO_HALO_
#endif
#define _2D_W_ _IRANGE_,_JRANGE_
#define _3D_W_ _IRANGE_,_JRANGE_,_KRANGE_

#define _2D_W_ _IRANGE_,_JRANGE_
#define _3D_W_ _IRANGE_,_JRANGE_,_KRANGE_

#define A1DFIELD  _KRANGE_
#define A2DFIELD  _IRANGE_HALO_,_JRANGE_HALO_
#define A3DFIELD  _IRANGE_HALO_,_JRANGE_HALO_,_KRANGE_

