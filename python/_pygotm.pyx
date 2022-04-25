# cython: language_level=3
# cython: profile=True

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

cdef extern void initialize(int nlev, const char* nml_file, double** ptke, double** ptkeo, double** peps, double** pL, double** pnum, double** pnuh) nogil
cdef extern void calculate(int nlev, double dt, double* h, double D, double taus, double taub, double z0s, double z0b, double* SS, double* NN) nogil
cdef extern void calculate_3d(int nx, int nz, int nz, int istart, int istop, int jstart, int jstop, double dt, int* mask, double* h3d, double* D, double* u_taus, double* u_taub, double* z0s, double* z0b, double* NN, double* SS, double* tke, double* tkeo, double* eps, double* L, double* num, double* nuh)
cdef extern void diff(int nlev, double dt, double cnpar, int posconc, double* h, int Bcup, int Bcdw, double Yup, double Ydw, double* nuY, double* Lsour, double* Qsour, double* Taur, double* Yobs, double* Y)
cdef extern void redirect_output(const char* stdout_file, const char* stderr_file) nogil
cdef extern void close_redirected_output() nogil

with tempfile.NamedTemporaryFile(delete=False) as f1, tempfile.NamedTemporaryFile(delete=False) as f2:
    stdout_path = f1.name
    stderr_path = f2.name
redirect_output(stdout_path.encode('ascii'), stderr_path.encode('ascii'))
stdout = open(stdout_path)
stderr = open(stderr_path)

def cleanup():
    close_redirected_output()
    stdout.close()
    stderr.close()
    os.remove(stdout_path)
    os.remove(stderr_path)
atexit.register(cleanup)

cdef class Mixing:
    cdef readonly numpy.ndarray tke, tkeo, eps, L, num, nuh
    cdef double* pnuh

    def __init__(self, int nlev, bytes nml_path=b''):
        cdef double* ptke
        cdef double* ptkeo
        cdef double* peps
        cdef double* pL
        cdef double* pnum
        cdef double* pnuh
        initialize(nlev, nml_path, &ptke, &ptkeo, &peps, &pL, &pnum, &pnuh)
        self.tke = numpy.asarray(<double[:nlev+1:1]> ptke)
        self.tkeo = numpy.asarray(<double[:nlev+1:1]> ptkeo)
        self.eps = numpy.asarray(<double[:nlev+1:1]> peps)
        self.L = numpy.asarray(<double[:nlev+1:1]> pL)
        self.num = numpy.asarray(<double[:nlev+1:1]> pnum)
        self.nuh = numpy.asarray(<double[:nlev+1:1]> pnuh)
        self.pnuh = pnuh

    def turbulence(self, double dt, const double[::1] h not None, double D, double u_taus, double u_taub, double z0s, double z0b, const double[::1] SS not None, const double[::1] NN not None):
        assert SS.shape[0] == NN.shape[0], 'Length of NN (%i) and SS (%i) should be identical.' % (NN.shape[0], SS.shape[0])
        assert h.shape[0] == NN.shape[0], 'Length of h (%i) should match that of NN and SS (%i).' % (h.shape[0], SS.shape[0])
        calculate(h.shape[0] - 1, dt, &h[0], D, u_taus, u_taub, z0s, z0b, &SS[0], &NN[0])

    def turbulence_3d(self, int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, double dt, const int[:, ::1] mask not None,
                      const double[:, :, ::1] h not None, const double[:, ::1] D not None, const double[:, ::1] u_taus not None, const double[:, ::1] u_taub not None,
                      const double[:, ::1] z0s not None, const double[:, ::1] z0b not None,  const double[:, :, ::1] NN not None,  const double[:, :, ::1] SS not None,
                      double[:, :, ::1] tke not None,  double[:, :, ::1] tkeo not None, double[:, :, ::1] eps not None,  double[:, :, ::1] L not None, double[:, :, ::1] num not None,  double[:, :, ::1] nuh not None):
        calculate_3d(nx, ny, nz, istart, istop, jstart, jstop, dt, &mask[0,0], &h[0,0,0], &D[0,0], &u_taus[0,0], &u_taub[0,0], &z0s[0,0], &z0b[0,0], &NN[0,0,0], &SS[0,0,0],
                     &tke[0,0,0], &tkeo[0,0,0], &eps[0,0,0], &L[0,0,0], &num[0,0,0], &nuh[0,0,0])

    def diffuse(self, double dt, const double[::1] h not None, double[::1] Y not None):
        cdef double[::1] Lsour, Qsour, Taur, Yobs
        assert h.shape[0] == Y.shape[0]
        assert h.shape[0] <= self.nuh.shape[0]
        Lsour = numpy.zeros_like(h)
        Qsour = numpy.zeros_like(h)
        Taur = numpy.full_like(h, 1e15)
        diff(h.shape[0] - 1, dt, 1., 0, &h[0], 1, 1, 0., 0., self.pnuh, &Lsour[0], &Qsour[0], &Taur[0], &Y[0], &Y[0])