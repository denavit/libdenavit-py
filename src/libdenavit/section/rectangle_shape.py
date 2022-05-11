from libdenavit.section import GeometricShape, FiberQuadPatch
from dataclasses import dataclass
from math import pi, tanh, sin, cos, sqrt
import numpy as np
import matplotlib.pyplot as plt
from libdenavit.design import available_strength

@dataclass
class Rectangle(GeometricShape):
    """
    Rectangle shape

    Parameters
    ----------
    H : float
        The width in the y direction
    B : float
        The width in the x direction
    rc : float
        Rectangle corner radius
    """

    H: float
    B: float
    rc: float = 0

    plot_fill_color = [0.90, 0.90, 0.90]

    @property
    def is_section_valid(self):
        tf1 = self.B > 0    # B should be positive
        tf2 = self.H > 0    # H should be positive
        tf3 = self.rc >= 0  # rc should be positive or zero
        tf4 = self.rc < min([self.B, self.H]) / 2
        tf = all([tf1, tf2, tf3, tf4])
        return tf

    def depth(self, axis):
        if axis.lower() == "x":
            d = self.H
        elif axis.lower() == "y":
            d = self.B
        else:
            raise ValueError(f'Unknown axis: {axis}')
        return d

    @property
    def A(self):
        a = self.H * self.B - (4 - pi) * self.rc ** 2
        return a

    @property
    def J(self):
        if self.rc != 0: raise ValueError("Not implemented for rounded rectangles")
        # Equation for J from Theory of Elasticity by Timoshenko and Goodier
        # (first two terms of the infinite series).
        # See also: Plaut, R. H., and Eatherton, M. R. (2017).
        # "Lateral-torsional buckling of butterfly-shaped beams with
        # rectangular cross section.” Engineering Structures, 136, 210–218.
        if self.H >= self.B:
            ar = self.H / self.B
            beta = 1.0 / 3.0 * (1 - 192.0 / pi ** 5 * 1 / ar * (tanh(pi * ar / 2.0) + tanh(3 * pi * ar / 2) / 243))
            j = beta * self.H * self.B ** 3
        else:
            ar = self.B / self.H
            beta = 1.0 / 3.0 * (1 - 192.0 / pi ** 5 * 1 / ar * (tanh(pi * ar / 2.0) + tanh(3 * pi * ar / 2) / 243))
            j = beta * self.B * self.H ** 3

        return j

    @property
    def Ix(self):
        if self.rc == 0:
            i = (1 / 12) * self.B * self.H ** 3
        else:
            i = ((1 / 12) * self.B * self.H ** 3
                 - 4 * ((1 / 12) * self.rc ** 4 + self.rc ** 2 * (self.H / 2 - self.rc / 2) ** 2)
                 + 4 * (pi / 16 - 4) * self.rc ** 4
                 + (pi / 4) * self.rc ** 2 * (self.H / 2 - (self.rc - (4 * self.rc) / (3 * pi))) ** 2)
        return i

    @property
    def Iy(self):
        if self.rc == 0:
            i = (1 / 12) * self.H * self.B ** 3
        else:
            i = ((1 / 12) * self.H * self.B ** 3
                 - 4 * ((1 / 12) * self.rc ** 4 + self.rc ** 2 * (self.B / 2 - self.rc / 2) ** 2)
                 + 4 * ((pi / 16 - 4 / (9 * pi)) * self.rc ** 4
                        + (pi / 4) * self.rc ** 2 * (self.B / 2 - (self.rc - (4 * self.rc) / (3 * pi))) ** 2))
        return i

    @property
    def Sx(self):
        return self.Ix/(0.5*self.H)

    @property
    def Sy(self):
        return self.Iy/(0.5*self.B)

    @property
    def Zx(self):
        if self.rc == 0:
            z = self.B * self.H ** 2 / 4
        else:
            z = (self.B * self.H ** 2 / 4
                 - 4 * ((1 - pi / 4) * self.rc ** 2) * (self.H / 2 - ((10 - 3 * pi) / (12 - 3 * pi)) * self.rc))
        return z

    @property
    def Zy(self):
        if self.rc == 0:
            z = self.H * self.B ** 2 / 4
        else:
            z = (self.H * self.B ** 2 / 4
                 - 4 * ((1 - pi / 4) * self.rc ** 2) * (self.B / 2 - ((10 - 3 * pi) / (12 - 3 * pi)) * self.rc))
        return z

    @property
    def boundary_points(self):
        if self.rc == 0:
            x = [self.B / 2, - self.B / 2, - self.B / 2, self.B / 2]
            y = [self.H / 2, self.H / 2, - self.H / 2, - self.H / 2]
            r = [0, 0, 0, 0]
        else:
            raise ValueError('Not yet implemented')

        return x, y, r

    def plot_section(self, line_width=2):
        if self.rc == 0:
            x = [self.B / 2, - self.B / 2, - self.B / 2, self.B / 2, self.B / 2]
            y = [self.H / 2, self.H / 2, - self.H / 2, - self.H / 2, self.H / 2]
        else:
            angles = np.linspace(0, pi / 2, 25)
            x = [(self.B / 2 + self.rc * cos(angles)),
                 (-self.B / 2 + self.rc * cos(angles + pi / 2)),
                 (-self.B / 2 + self.rc * cos(angles + pi)),
                 (self.B / 2 + self.rc * cos(angles + 1.5 * pi)), self.B / 2 + self.rc]
            y = [(self.H / 2 + self.rc * sin(angles)),
                 (self.H / 2 + self.rc * sin(angles + pi / 2)),
                 (-self.H / 2 + self.rc * sin(angles + pi)),
                 (-self.H / 2 + self.rc * sin(angles + 1.5 * pi)), self.H / 2]

        plt.plot(x, y)

    def add_to_fiber_section(self, fiber_section, mat_id):
        if self.rc == 0:
            a = FiberQuadPatch(-self.B / 2, - self.H / 2,
                               -self.B / 2,   self.H / 2,
                                self.B / 2,   self.H / 2,
                                self.B / 2,  -self.H / 2, mat_id)
            fiber_section.add_fibers(a)
        else:
            ValueError('Not yet implemented')


class PlateMember_AISC2016:
    def __init__(self,section,Fy,E,Lcx,Lcy,strength_type):
        self.section = section
        self.Fy = Fy
        self.E = E
        self.Lcx = Lcx
        self.Lcy = Lcy
        self.strength_type = strength_type

    def Pnt(self):
        Pn = self.Fy*self.section.A
        return available_strength(Pn,self.strength_type,0.9,1.67)

    def Pnc(self):
        Lc_over_r = max(self.Lcx/self.section.rx,self.Lcy/self.section.ry)
        
        if Lc_over_r == 0.0:
            Fcr = self.Fy
        else:
            Fe = pi**2*self.E/Lc_over_r**2 
            if Lc_over_r <= 4.71*sqrt(self.E/self.Fy):
                Fcr = (0.658**(self.Fy/Fe))*self.Fy
            else:
                Fcr = 0.877*Fe
                
        Pn = Fcr*self.section.A
        return available_strength(Pn,self.strength_type,0.9,1.67)