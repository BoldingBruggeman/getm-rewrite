# cython: language_level=3
# cython: profile=True

cimport cython

cimport numpy
import numpy

cdef extern void initialize(int nlev, const char* inputdir, double** ptke, double** ptkeo, double** peps, double** pL, double** pnum, double** pnuh) nogil
cdef extern void calculate(int nlev, double dt, double* h, double D, double taus, double taub, double z0s, double z0b, double* SS, double* NN) nogil
cdef extern void diff(int nlev, double dt, double cnpar, int posconc, double* h, int Bcup, int Bcdw, double Yup, double Ydw, double* nuY, double* Lsour, double* Qsour, double* Taur, double* Yobs, double* Y)

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

    def diffuse(self, double dt, const double[::1] h not None, double[::1] Y not None):
        cdef double[::1] Lsour, Qsour, Taur, Yobs
        assert h.shape[0] == Y.shape[0]
        assert h.shape[0] <= self.nuh.shape[0]
        Lsour = numpy.zeros_like(h)
        Qsour = numpy.zeros_like(h)
        Taur = numpy.full_like(h, 1e15)
        diff(h.shape[0] - 1, dt, 1., 0, &h[0], 1, 1, 0., 0., self.pnuh, &Lsour[0], &Qsour[0], &Taur[0], &Y[0], &Y[0])