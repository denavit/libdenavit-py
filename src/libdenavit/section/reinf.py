from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt
from libdenavit.section import FiberSingle


class Reinf:
    Ab = None
    _db = None
    color_steelFill = [0.75, 0.75, 0.75]

    def I(self, axis):
        x, y = self.coordinates
        if axis.lower() == "x":
            i = sum(self.Ab * y ** 2)
        elif axis.lower == "y":
            i = sum(self.Ab * x ** 2)
        else:
            raise ValueError("Unknown Axis")

        return i

    def plot_section(self, line_width=2):
        xb, yb = self.coordinates
        plt.plot(xb, yb, 'o', color=self.color_steelFill, linewidth=line_width, )


@dataclass
class ReinfRect(Reinf):
    """
    A ReinfRect object is a rebar layer for a rectangular section with a given width and height.

    Parameters
    ----------
    Bx : float
        The width of the rebar layer in x direction.
    By : float
        The height of the rebar layer in y direction.
    nbx : int
        number of bars in the x direction of the section.
    nby : int
        number of bars in the y direction of the section.
    Ab : float
        The area of the rebar layer.
    xc: float
        The x coordinate of the center of the rebar layer.
    yc: float
        The y coordinate of the center of the rebar layer.
    """
    Bx: float
    By: float
    nbx: int
    nby: int
    Ab: float
    xc: float = 0
    yc: float = 0

    @property
    def coordinates(self):
        x1 = np.linspace(-0.5 * self.Bx + self.xc, +0.5 * self.Bx + self.xc, self.nbx)
        x2 = np.linspace(+0.5 * self.Bx + self.xc, +0.5 * self.Bx + self.xc, self.nby)
        x3 = np.linspace(+0.5 * self.Bx + self.xc, -0.5 * self.Bx + self.xc, self.nbx)
        x4 = np.linspace(-0.5 * self.Bx + self.xc, -0.5 * self.Bx + self.xc, self.nby)
        y1 = np.linspace(+0.5 * self.By + self.yc, +0.5 * self.By + self.yc, self.nbx)
        y2 = np.linspace(+0.5 * self.By + self.yc, -0.5 * self.By + self.yc, self.nby)
        y3 = np.linspace(-0.5 * self.By + self.yc, -0.5 * self.By + self.yc, self.nbx)
        y4 = np.linspace(-0.5 * self.By + self.yc, +0.5 * self.By + self.yc, self.nby)
        x = np.concatenate((x1[0:-1], x2[0:-1], x3[0:-1], x4[0:-1]))
        y = np.concatenate((y1[0:-1], y2[0:-1], y3[0:-1], y4[0:-1]))
        return x, y

    @property
    def boundaries(self):
        x = [-0.5 * self.Bx + self.xc,
             -0.5 * self.Bx + self.xc,
             +0.5 * self.Bx - self.xc,
             +0.5 * self.Bx - self.xc, ]

        y = [+0.5 * self.By - self.yc,
             -0.5 * self.By + self.yc,
             +0.5 * self.By - self.yc,
             -0.5 * self.By + self.yc, ]
        return x, y

    @property
    def num_bars(self):
        n = 2 * (self.nbx + self.nby) - 4
        return n

    @property
    def db(self):
        if self._db is not None:
            return self._db
        else:
            return sqrt(4 * self.Ab / pi)

    @db.setter
    def db(self, x):
        self._db = x

    def add_to_fiber_section(self, fiber_section, id_reinf, id_conc):
        [x, y] = self.coordinates
        for j in range(len(x)):
            a = FiberSingle(self.Ab, x[j], y[j], id_reinf, id_conc)
            fiber_section.add_fibers(a)
