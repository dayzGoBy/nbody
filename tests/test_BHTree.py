import numpy as np
from ..barnes_hut.BHTree import Tracker, Body


def cmp(a, b, eps=0.001):
    return abs(a - b) / abs(b) < eps


def test_solar():
    SIZE = np.float64(2.5e12)
    BODIES = [Body(1.4960e+11 + 0j,  0 + 2.9800e+04j,  5.9740e+24),
              Body(2.2790e+11 + 0j, 0 + 2.4100e+04j, 6.4190e+23),
              Body(5.7900e+10 + 0j, 0 + 4.7900e+04j, 3.3020e+23),
              Body(0, 0, 1.9890e+30),
              Body(1.0820e+11, 3.5000e+04j, 4.8690e+24)]

    T = Tracker(BODIES, SIZE)

    T.update(N=1000)

    assert cmp(T.get(0), 1.496e+11 + 2.98e+07j)
    assert cmp(T.get(3), 0.0264958 + 3.00701e-06j)
    assert cmp(T.get(4), 1.082e+11 + 3.5e+07j)


def test_big():
    SIZE = np.float64(2.5e12)
    BODIES = [Body(0, 0, 2e32),
              Body(6.25e10, 5.66e05j, 6e24),
              Body(-6.25e10, -5.66e05j, 6e24),
              Body(6.250e10j, -5.66e05, 6e24),
              Body(-6.250e10j, 5.66e05, 6e24),
              Body(2.828e11 + 2.828e11j, 100000 - 100000j, 1e26),
              Body(2.828e11 - 2.828e11j, -100000 - 100000j,  1e26),
              Body(-2.828e11 - 2.828e11j, -100000 + 100000j,   1e26),
              Body(-2.828e11 + 2.828e11j, 100000 + 100000j,   1e26),
              Body(-1.828e11 + 1.828e11j, 100000 + 100000j,   1e16)]

    T = Tracker(BODIES.copy(), SIZE)

    T.update(N=1000)

    assert cmp(T.get(2), -6.24983e+10 - 5.65995e+08j)
    assert cmp(T.get(4), 5.65995e+08 - 6.24983e+10j)
