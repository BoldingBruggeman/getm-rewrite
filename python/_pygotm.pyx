# cython: language_level=3
# cython: profile=True

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

cdef extern void initialize(int nlev, const char* nml_file, const char* yaml_file, double** ptke, double** peps, double** pL, double** pnum, double** pnuh) nogil
cdef extern void finalize() nogil
cdef extern void calculate(int nlev, double dt, const double* h, double D, double taus, double taub, double z0s, double z0b, const double* SS, const double* NN) nogil
cdef extern void calculate_3d(int nx, int nz, int nz, int istart, int istop, int jstart, int jstop, double dt, const int* mask, const double* h3d, const double* D, const double* u_taus, const double* u_taub, const double* z0s, const double* z0b, const double* NN, const double* SS, double* tke, double* eps, double* L, double* num, double* nuh)
cdef extern void diff(int nlev, double dt, double cnpar, int posconc, const double* h, int Bcup, int Bcdw, double Yup, double Ydw, const double* nuY, const double* Lsour, const double* Qsour, const double* Taur, const double* Yobs, double* Y)
cdef extern void redirect_output(const char* stdout_file, const char* stderr_file) nogil
cdef extern void close_redirected_output() nogil

with tempfile.NamedTemporaryFile(delete=False) as f1, tempfile.NamedTemporaryFile(delete=False) as f2:
    stdout_path = f1.name
    stderr_path = f2.name
redirect_output(stdout_path.encode('ascii'), stderr_path.encode('ascii'))
stdout = open(stdout_path)
stderr = open(stderr_path)

def cleanup():
    finalize()
    close_redirected_output()
    stdout.close()
    stderr.close()
    os.remove(stdout_path)
    os.remove(stderr_path)
atexit.register(cleanup)

cdef class Mixing:
    cdef readonly numpy.ndarray tke, eps, L, num, nuh
    cdef double* pnuh

    def __init__(self, int nlev, bytes nml_path=b'', bytes yaml_path=b''):
        cdef double* ptke
        cdef double* peps
        cdef double* pL
        cdef double* pnum
        cdef double* pnuh
        initialize(nlev, nml_path, yaml_path, &ptke, &peps, &pL, &pnum, &pnuh)
        self.tke = numpy.asarray(<double[:nlev+1:1]> ptke)
        self.eps = numpy.asarray(<double[:nlev+1:1]> peps)
        self.L = numpy.asarray(<double[:nlev+1:1]> pL)
        self.num = numpy.asarray(<double[:nlev+1:1]> pnum)
        self.nuh = numpy.asarray(<double[:nlev+1:1]> pnuh)
        self.pnuh = pnuh

    def turbulence(self, double dt, const double[::1] h not None, double D, double u_taus, double u_taub, double z0s, double z0b, const double[::1] SS not None, const double[::1] NN not None):
        assert SS.shape[0] == NN.shape[0], 'Length of NN (%i) and SS (%i) should be identical.' % (NN.shape[0], SS.shape[0])
        assert h.shape[0] == NN.shape[0], 'Length of h (%i) should match that of NN and SS (%i).' % (h.shape[0], SS.shape[0])
        calculate(<int>h.shape[0] - 1, dt, &h[0], D, u_taus, u_taub, z0s, z0b, &SS[0], &NN[0])

    def turbulence_3d(self, int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, double dt, const int[:, ::1] mask not None,
                      const double[:, :, ::1] h not None, const double[:, ::1] D not None, const double[:, ::1] u_taus not None, const double[:, ::1] u_taub not None,
                      const double[:, ::1] z0s not None, const double[:, ::1] z0b not None,  const double[:, :, ::1] NN not None,  const double[:, :, ::1] SS not None,
                      double[:, :, ::1] tke not None,  double[:, :, ::1] eps not None,  double[:, :, ::1] L not None, double[:, :, ::1] num not None,  double[:, :, ::1] nuh not None):
        assert mask.shape[1] == nx and mask.shape[0] == ny
        assert h.shape[2] == nx and h.shape[1] == ny and h.shape[0] == nz
        assert D.shape[1] == nx and D.shape[0] == ny
        assert u_taus.shape[1] == nx and u_taus.shape[0] == ny
        assert u_taub.shape[1] == nx and u_taub.shape[0] == ny
        assert z0s.shape[1] == nx and z0s.shape[0] == ny
        assert z0b.shape[1] == nx and z0b.shape[0] == ny
        assert NN.shape[2] == nx and NN.shape[1] == ny and NN.shape[0] == nz + 1
        assert SS.shape[2] == nx and SS.shape[1] == ny and SS.shape[0] == nz + 1
        assert tke.shape[2] == nx and tke.shape[1] == ny and tke.shape[0] == nz + 1
        assert eps.shape[2] == nx and eps.shape[1] == ny and eps.shape[0] == nz + 1
        assert L.shape[2] == nx and L.shape[1] == ny and L.shape[0] == nz + 1
        assert num.shape[2] == nx and num.shape[1] == ny and num.shape[0] == nz + 1
        assert nuh.shape[2] == nx and nuh.shape[1] == ny and nuh.shape[0] == nz + 1

        # note we convert from 0-based start indices (Python/C) to 1-based start indices (Fortran)
        # stop indices are already fine because Python uses the first index that is excluded, whereas Fortran uses the last index that is included. These happen to be the same.
        calculate_3d(nx, ny, nz, istart + 1, istop, jstart + 1, jstop, dt, &mask[0,0], &h[0,0,0], &D[0,0], &u_taus[0,0], &u_taub[0,0], &z0s[0,0], &z0b[0,0], &NN[0,0,0], &SS[0,0,0],
                     &tke[0,0,0], &eps[0,0,0], &L[0,0,0], &num[0,0,0], &nuh[0,0,0])

    def diffuse(self, double dt, const double[::1] h not None, double[::1] Y not None):
        cdef double[::1] Lsour, Qsour, Taur, Yobs
        assert h.shape[0] == Y.shape[0]
        assert h.shape[0] <= self.nuh.shape[0]
        Lsour = numpy.zeros_like(h)
        Qsour = numpy.zeros_like(h)
        Taur = numpy.full_like(h, 1e15)
        diff(<int>h.shape[0] - 1, dt, 1., 0, &h[0], 1, 1, 0., 0., self.pnuh, &Lsour[0], &Qsour[0], &Taur[0], &Y[0], &Y[0])