# cython: language_level=3

cimport cython

cimport numpy
import numpy

cdef extern void initialize(int nlev, const char* inputdir, double** ptke, double** peps, double** pL, double** pnum, double** pnuh) nogil
cdef extern void calculate(int nlev, double dt, double* h, double D, double taus, double taub, double z0s, double z0b, double* SS, double* NN) nogil
cdef extern void diff(int nlev, double dt, double cnpar, int posconc, double* h, int Bcup, int Bcdw, double Yup, double Ydw, double* nuY, double* Lsour, double* Qsour, double* Taur, double* Yobs, double* Y)

cdef class Mixing:
    cdef readonly numpy.ndarray tke, eps, L, num, nuh
    cdef double* pnuh

    def __init__(self, int nlev, bytes inputdir=b'.'):
        cdef double* ptke
        cdef double* peps
        cdef double* pL
        cdef double* pnum
        cdef double* pnuh
        initialize(nlev, inputdir, &ptke, &peps, &pL, &pnum, &pnuh)
        self.tke = numpy.asarray(<double[:nlev+1:1]> ptke)
        self.eps = numpy.asarray(<double[:nlev+1:1]> peps)
        self.L = numpy.asarray(<double[:nlev+1:1]> pL)
        self.num = numpy.asarray(<double[:nlev+1:1]> pnum)
        self.nuh = numpy.asarray(<double[:nlev+1:1]> pnuh)
        self.pnuh = pnuh

    def turbulence(self, int nlev, double dt, const double[::1] h not None, double D, double taus, double taub, double z0s, double z0b, const double[::1] SS not None, const double[::1] NN not None):
        calculate(nlev, dt, &h[0], D, taus, taub, z0s, z0b, &SS[0], &NN[0])

    def diffuse(self, int nlev, double dt, const double[::1] h not None, double[::1] Y not None):
        cdef double[::1] Lsour, Qsour, Taur, Yobs
        Lsour = numpy.zeros_like(h)
        Qsour = numpy.zeros_like(h)
        Taur = numpy.full_like(h, 1e15)
        diff(nlev, dt, 1., 0, &h[0], 1, 1, 0., 0., self.pnuh, &Lsour[0], &Qsour[0], &Taur[0], &Y[0], &Y[0])