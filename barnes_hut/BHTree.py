import numpy as np
from enum import IntEnum
import abc
from collections import deque


class DIR(IntEnum):
    NULL = -1
    SE = 0
    SW = 1
    NE = 2
    NW = 3


class Body:
    def __init__(self, position, velocity, mass):
        self.__mass = mass * np.float64(6.674184e-11)
        self.__pos = position
        self.__vel = velocity
        self.__acc = np.complex128(0)

    def get_mass(self):
        return self.__mass

    def get_pos(self):
        return self.__pos

    def get_vel(self):
        return self.__pos

    def add_force(self, pos, mass):
        r = pos - self.__pos
        len_r = abs(r)
        self.__acc += r * mass / (len_r * len_r * len_r)

    def add_force_body(self, other):
        self.add_force(other.__pos, other.__mass)

    def update(self, dt):
        dv = self.__acc * dt
        self.__pos += self.__vel * dt + dv * dt / 2
        self.__vel += dv
        self.__acc = np.complex128(0)


class Quad:
    def __init__(self, center, length):
        self.center = center
        self.length = length
        self.__half_len = length / 2

    def contains(self, p):
        return (p.real < (self.center.real + self.__half_len)
                and p.real > (self.center.real - self.__half_len)
                and p.imag < (self.center.imag + self.__half_len)
                and p.imag > (self.center.imag - self.__half_len))

    def subquad(self, dir: DIR):
        if dir == -1:
            raise ValueError("incorrect direction")
        qr = self.__half_len / 2
        mul_x = -1 if dir == DIR.SW or dir == DIR.NW else 1
        mul_y = 1 if dir == DIR.NE or dir == DIR.NW else -1
        return Quad(self.center +
                    np.complex128(qr * mul_x + 1j * qr * mul_y),
                    self.__half_len)


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

    def contains(self, body):
        return self.quad.contains(body.get_pos())


class LNode(Node):
    def __init__(self, quad, parent, dir, body):
        Node.__init__(self, quad, parent, dir)
        self.body = body

    @classmethod
    def create_from_parent(cls, parent, dir, body=None):
        return cls(parent.quad.subquad(dir), parent, dir, body)

    @classmethod
    def create_from_pos(cls, quad, body):
        return cls(quad, None, None, body)

    def is_leaf(self):
        return True


class INode(Node):
    def __init__(self, quad, parent, dir):
        Node.__init__(self, quad, parent, dir)
        self.children = [LNode.create_from_parent(self, DIR.SE),
                         LNode.create_from_parent(self, DIR.SW),
                         LNode.create_from_parent(self, DIR.NE),
                         LNode.create_from_parent(self, DIR.NW)]
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
        self.children[dir] = INode.create_from_parent(self, dir)

    def where(self, pos):
        is_w = pos.real < self.quad.center.real
        is_n = pos.imag > self.quad.center.imag
        return self.children[is_n * 2 + is_w]

    def get(self, dir):
        return self.children[dir]

    def calc_mass(self):
        self.mass = np.float64(0)
        for ch in self.children:
            if ch.is_leaf() and ch.body is not None:
                self.mass += 0 if ch.body is None else ch.body.get_mass()
            elif not ch.is_leaf():
                self.mass += ch.calc_mass()
        return self.mass

    def calc_mass_center(self):
        self.mass_center = np.complex128(0)
        for ch in self.children:
            if ch.is_leaf() and ch.body is not None:
                self.mass_center += (ch.body.get_mass() * ch.body.get_pos() /
                                     self.mass)
            elif not ch.is_leaf():
                ch.calc_mass_center()
                self.mass_center += ch.mass * ch.mass_center / self.mass
        return self.mass_center


class BHTree:
    THETA = 0.5

    def __init__(self, bodies, size):
        self.root = None
        self.size = size
        for b in bodies:
            self.insert(b)
        if len(bodies) >= 2:
            self.root.calc_mass()
            self.root.calc_mass_center()

    def insert(self, body):
        if self.root is None:
            self.root = LNode.create_from_pos(
                Quad(0, self.size), body)
            return

        ptr = self.root
        while not ptr.is_leaf():
            ptr = ptr.where(body.get_pos())

        if ptr.body is None:
            ptr.body = body
            return

        another = ptr.body
        dir = ptr.dir
        split_node = None
        parent = ptr.parent

        if ptr.parent is None:
            self.root = INode.create_from_pos(ptr.quad)
            split_node = self.root
        else:
            parent.set_internal(dir)
            split_node = parent

        self.split(body, another, split_node)

    def split(self, one, two, cur):
        one_p = one.get_pos()
        two_p = two.get_pos()
        one_l = cur.where(one_p)
        two_l = cur.where(two_p)

        if one_l == two_l:
            dir = one_l.dir
            cur.set_internal(dir)
            self.split(one, two, cur.get(dir))
        else:
            one_l.body = one
            two_l.body = two

    def calc_force(self, body):
        q = deque()
        q.append(self.root)

        while len(q) > 0:
            node = q.popleft()

            if (node.is_leaf() and node.body is not None and
                    node.body is not body):
                body.add_force_body(node.body)
            elif not node.is_leaf() and node.contains(body):
                for ch in node.children:
                    q.append(ch)
            elif not node.is_leaf() and not node.contains(body):
                theta = node.quad.length / abs(body.get_pos()
                                               - node.mass_center)
                if theta < self.THETA:
                    body.add_force(node.mass_center, node.mass)
                else:
                    for ch in node.children:
                        q.append(ch)


class Tracker:
    def __init__(self, bodies, size, dt=np.float64(1)):
        self.__bodies = bodies.copy()
        self.__size = size
        self.dt = dt

    def update(self, N=1):
        for _ in range(N):
            tree = BHTree(self.__bodies, self.__size)
            for b in self.__bodies:
                tree.calc_force(b)

            for b in self.__bodies:
                b.update(self.dt)

    def get(self, index):
        return self.__bodies[index].get_pos()

    def get_all(self):
        return np.array(list(map(lambda x: x.get_pos(), self.__bodies)))

    def get_all_pairs(self):
        return np.array(list(map(
            lambda x: (x.get_pos().real, x.get_pos().imag,),
            self.__bodies)))
