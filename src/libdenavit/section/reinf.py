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
        elif axis.lower() == "y":
            i = sum(self.Ab * x ** 2)
        else:
            raise ValueError("Unknown Axis")

        return i

    def plot_section(self, **kwargs):
        xb, yb = self.coordinates
        plt.plot(xb, yb, 'o', **kwargs)
    
    def add_to_fiber_section(self, fiber_section, id_reinf, id_conc):
        [x, y] = self.coordinates
        for j in range(len(x)):
            a = FiberSingle(self.Ab, x[j], y[j], id_reinf, id_conc)
            fiber_section.add_fibers(a)

    @property
    def db(self):
        if self._db is not None:
            return self._db
        else:
            return sqrt(4 * self.Ab / pi)

    @db.setter
    def db(self, x):
        self._db = x


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
    def num_bars(self):
        n = 2 * (self.nbx + self.nby) - 4
        return n


@dataclass
class ReinfCirc(Reinf):
    """
    A ReinfCirc object is a rebar layer for a circular section with a given radius.

    Parameters
    ----------
    rc : float
        radius from center of the bar pattern to the center of the bar.
    num_bars : int
        number of bars in the circumference.
    Ab : float
        area of each rebar.
    xc : float
        The x coordinate of the center of the rebar layer.
    yc : float
        The y coordinate of the center of the rebar layer.

    """
    rc: float
    num_bars: int
    Ab: float
    xc: float = 0
    yc: float = 0

    @property
    def coordinates(self):
        angles = np.linspace(0, 2 * np.pi, self.num_bars + 1)[:-1]
        x = self.xc + self.rc * np.cos(angles)
        y = self.yc + self.rc * np.sin(angles)
        return x, y


class ReinfIntersectingLoops(Reinf):
    """
    A ReinfIntersectingLoops object is a rebar layer for a obround section with two circular patterns.

    Parameters
    ----------
    D  : float
        diameter of the intersecting loops (measured to the center of the bars).
    a  : float
        seperation of the centers of the two circlar bar patterns.
    num_bars : int
        number of bars in the patter, must be even.
    Ab : float
        area of each reinforcing bar.
    xc : float
        The x coordinate of the center of the whole pattern.
    yc : float
        The y coordinate of the center of the whole pattern.

    """

    def __init__(self,D,a,num_bars,Ab,xc=0,yc=0):
        self.D = D
        self.a = a
        self.num_bars = num_bars
        self.Ab = Ab
        self.xc = xc
        self.yc = yc

    @property
    def coordinates(self):
        if self.a >= self.D:
            raise ValueError(f'a should be less than D')
        if self.num_bars % 2 != 0:
            raise ValueError(f'num_bars should be even ({self.num_bars = })')
        
        theta = np.arccos(self.a/self.D)
        
        angles = np.linspace(-np.pi+theta, np.pi-theta, self.num_bars//2+1)
        x1 = self.xc + self.D/2 * np.cos(angles) + self.a/2
        y1 = self.yc + self.D/2 * np.sin(angles)
        
        angles = np.linspace(theta, 2*np.pi-theta, self.num_bars//2+1)[1:-1]
        x2 = self.xc + self.D/2 * np.cos(angles) - self.a/2
        y2 = self.yc + self.D/2 * np.sin(angles)
        
        x = np.append(x1,x2)
        y = np.append(y1,y2)
        
        return x, y