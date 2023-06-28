import numpy as np
from numba import njit
from enum import IntEnum
import abc


class DIR(IntEnum):
    NULL = -1
    SE = 0
    SW = 1
    NE = 2
    NW = 3


class Body:
    G = np.float64(6.674e-11)

    def __init__(self, mass: np.float64, position: np.complex128,
                 velocity: np.complex128):
        self.__mass = mass / self.G
        self.__vel = velocity
        self.__pos = position
        self.__acc = np.complex128(0)

    def get_mass(self):
        return self.__mass

    def get_pos(self):
        return self.__pos

    @njit(inline='always')
    def add_force(self, other):
        r = other.pos - self.__pos
        len_r = abs(r)
        self.__acc += r * other.__mass / (len_r * len_r * len_r)

    @njit(inline='always')
    def update(self, dt):
        dv = self.acc * dt
        self.__pos += self.__vel * dt + dv * dt / 2
        self.__vel += dv
        self.acc = np.complex128(0)


class Quad:
    def __init__(self, center: np.complex128, lenght: np.float64):
        self.center = center
        self.length = lenght
        self.__half_len = self.length / 2

    def contains(self, p: np.complex128):
        return (p.real <= self.center.real + self.__half_len
                and p.real >= self.center.real - self.__half_len
                and p.imag <= self.center.imag + self.__half_len
                and p.imag >= self.center.imag - self.__half_len)

    def subquad(self, dir: DIR):
        if dir == -1:
            raise ValueError("incorrect direction")
        qr = self.__half_len / 2
        mul_x = -1 if dir == DIR.SW or dir == DIR.NW else 1
        mul_y = 1 if dir == DIR.NE or dir == DIR.NW else -1
        return Quad(self.center + np.complex128(qr * mul_x + 1j * qr * mul_y),
                    self.__half_len)


class BHTree:
    class Node(metaclass=abc.ABCMeta):
        def __init__(self, quad, parent, dir):
            self.quad = quad
            self.parent = parent
            self.dir = dir

        @classmethod
        def create_from_parent(cls, parent, dir):
            return cls(parent.quad.subquad(dir), parent, dir)

        @classmethod
        def create_from_pos(cls, quad):
            return cls(quad, None, None)

        @abc.abstractmethod
        def is_leaf(self):
            pass

    class INode(Node):
        def __init__(self, quad, parent, dir):
            super.__init__(quad, parent, dir)
            self.children = [None] * 4
            self.mass = np.float64(0)
            self.mass_center = np.complex128(0)

        @classmethod
        def create_from_parent(cls, parent, dir):
            return cls(parent.quad.subquad(dir), parent, dir)

        @classmethod
        def create_from_pos(cls, quad):
            return cls(quad, None, None)

        def is_leaf(self):
            return False

        def set_internal(self, dir):
            # self.children[dir] = INode()
            pass

        def where(self, pos):
            is_w = pos.real <= self.quad.center.real
            is_n = pos.imag >= self.quad.center.imag
            return self.children[is_n * 2 + is_w]

        def get(self, dir):
            return self.children[dir]

        def calc_mass(self):
            self.mass = np.float64(0)
            for ch in self.children:
                if ch.is_leaf():
                    self.mass += 0 if ch.body is None else ch.mass
                else:
                    self.mass += ch.calc_mass()
            return self.mass

        def calc_mass_center(self):
            self.mass_center = np.complex128(0)
            for ch in self.children:
                if ch.is_leaf():
                    if ch.body is not None:
                        self.mass_center += (ch.body.mass *
                                             ch.body.pos / self.mass)
                else:
                    self.mass += ch.calc_mass_center()
            return self.mass_center

    class LNode(Node):
        def __init__(self, quad, parent, dir, body):
            super.__init__(quad, parent, dir)
            self.body = body

        @classmethod
        def create_from_parent(cls, parent, dir, body=None):
            return cls(parent.quad.subquad(dir), parent, dir, body)

        @classmethod
        def create_from_pos(cls, quad, body):
            return cls(quad, None, None, body)

        def is_leaf(self):
            return True

    def __init__(self, bodies):
        self.root = None
        for b in bodies:
            self.insert(b)
        if len(bodies) >= 2:
            # TODO: calc mass and mass center
            pass

    def insert(self, body):
        pass

    def calc_force(self, body):
        pass


class Tracker:
    def __init__(self, bodies, dt=np.float64(0.1)):
        self.__bodies = bodies.copy()
        self.dt = dt

    @njit
    def update(self, N=1):
        tree = BHTree(self.__bodies)

        for i in range(N):
            for b in self.__bodies:
                tree.calc_force(b)

        for b in self.__bodies:
            b.update()
