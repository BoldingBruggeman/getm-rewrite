# cython: language_level=3
# cython: profile=True

cimport cython

cimport numpy
import numpy

cdef extern void solar_zenith_angle_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, int yday, double hh, const double* dlon,  const double* dlat,  double* zenith_angle)
cdef extern void shortwave_radiation_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, int yday, const double* zenith_angle, const double* dlat, const double* cloud, double* swr)
cdef extern void longwave_radiation_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, int method, const double* dlat, const double* tw, const double* ta, const double* cloud, const double* ea, const double* qa, double* ql)
cdef extern void humidity_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, int method, const double* hum, const double* airp, const double* tw, const double* ta, double* es, double* ea, double* qs, double* qa, double* rhoa)
cdef extern void albedo_water_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, int method, const double* zen, int yday, double* albedo)
cdef extern void kondo_2d(int nx, int ny, int imin, int imax, int jmin, int jmax, const double* sst, const double* airt, const double* w, double* cd_mom, double* cd_sensible, double* cd_latent)

def solar_zenith_angle(int yday, double hh, const double[:, ::1] dlon,  const double[:, ::1] dlat, double[:, ::1] zenith_angle, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = dlon.shape[0]
    if jmax is None:
        jmax = dlon.shape[1]
    solar_zenith_angle_2d(<int>dlon.shape[1], <int>dlon.shape[0], jmin + 1, jmax, imin + 1, imax, yday, hh, &dlon[0,0], &dlat[0,0], &zenith_angle[0,0])

def shortwave_radiation(int yday, const double[:, ::1] zenith_angle, const double[:, ::1] dlat, const double[:, ::1] cloud, double[:, ::1] swr, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = dlat.shape[0]
    if jmax is None:
        jmax = dlat.shape[1]
    shortwave_radiation_2d(<int>dlat.shape[1], <int>dlat.shape[0], jmin + 1, jmax, imin + 1, imax, yday, &zenith_angle[0,0], &dlat[0,0], &cloud[0,0], &swr[0,0])

def humidity(int method, const double[:, ::1] hum, const double[:, ::1] airp,  const double[:, ::1] tw, const double[:, ::1] ta, double[:, ::1] es, double[:, ::1] ea, double[:, ::1] qs, double[:, ::1] qa, double[:, ::1] rhoa, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = hum.shape[0]
    if jmax is None:
        jmax = hum.shape[1]
    humidity_2d(<int>hum.shape[1], <int>hum.shape[0], jmin + 1, jmax, imin + 1, imax, method, &hum[0,0], &airp[0,0], &tw[0,0], &ta[0,0], &es[0,0], &ea[0,0], &qs[0,0], &qa[0,0], &rhoa[0,0])

def albedo_water(int method, const double[:, ::1] zen, int yday, double[:, ::1] albedo, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = zen.shape[0]
    if jmax is None:
        jmax = zen.shape[1]
    albedo_water_2d(<int>zen.shape[1], <int>zen.shape[0], jmin + 1, jmax, imin + 1, imax, method, &zen[0,0], yday, &albedo[0,0])

def longwave_radiation(int method, const double[:, ::1] dlat, const double[:, ::1] tw, const double[:, ::1] ta, const double[:, ::1] cloud, const double[:, ::1] ea, const double[:, ::1] qa, double[:, ::1] ql, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = dlat.shape[0]
    if jmax is None:
        jmax = dlat.shape[1]
    assert tw.shape[0] == dlat.shape[0] and tw.shape[1] == dlat.shape[1]
    assert ta.shape[0] == dlat.shape[0] and ta.shape[1] == dlat.shape[1]
    assert cloud.shape[0] == dlat.shape[0] and cloud.shape[1] == dlat.shape[1]
    assert ea.shape[0] == dlat.shape[0] and ea.shape[1] == dlat.shape[1]
    assert qa.shape[0] == dlat.shape[0] and qa.shape[1] == dlat.shape[1]
    assert ql.shape[0] == dlat.shape[0] and ql.shape[1] == dlat.shape[1]
    longwave_radiation_2d(<int>dlat.shape[1], <int>dlat.shape[0], jmin + 1, jmax, imin + 1, imax, method, &dlat[0,0], &tw[0,0], &ta[0,0], &cloud[0,0], &ea[0,0], &qa[0,0], &ql[0,0])

def transfer_coefficients(int method, const double[:, ::1] tw, const double[:, ::1] ta, const double[:, ::1] w, double[:, ::1] cd_mom, double[:, ::1] cd_sensible, double[:, ::1] cd_latent, int imin=0, imax=None, int jmin=0, jmax=None):
    if imax is None:
        imax = tw.shape[0]
    if jmax is None:
        jmax = tw.shape[1]
    assert tw.shape[0] == ta.shape[0] and tw.shape[1] == ta.shape[1]
    assert tw.shape[0] == w.shape[0] and tw.shape[1] == w.shape[1]
    assert tw.shape[0] == cd_mom.shape[0] and tw.shape[1] == cd_mom.shape[1]
    assert tw.shape[0] == cd_sensible.shape[0] and tw.shape[1] == cd_sensible.shape[1]
    assert tw.shape[0] == cd_latent.shape[0] and tw.shape[1] == cd_latent.shape[1]
    kondo_2d(<int>tw.shape[1], <int>tw.shape[0], jmin + 1, jmax, imin + 1, imax, &tw[0,0], &ta[0,0], &w[0,0], &cd_mom[0,0], &cd_sensible[0,0], &cd_latent[0,0])
