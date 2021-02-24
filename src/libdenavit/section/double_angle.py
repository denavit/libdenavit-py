import dataclasses
from math import sqrt,atan,radians,tan

from . import database


# Dataclasses give some nice benefits for simple classes, like an automatic
# __repr__. So when printing the object in the console you'll get 
#
#   DoubleAngle(d=..., b=..., t=..., s=..., name=...)
#
# instead of
#
#   "<pystructe.section.double_angle.DoubleAngle at [address]>"
#
# See https://docs.python.org/3/library/dataclasses.html for more. Requires
# Python 3.7+.
#
@dataclasses.dataclass
class DoubleAngle:
    """
    Parameters
    ----------
    d : float
        leg length along the y-axis
    b : float
        leg length along the x-axis
    t : float
        thickness
    s : float
        spacing between angles
    name : str, optional
        name of the section

    Note that the definition of d and b are swapped from the Angle class.

    Calculations neglect the toe fillet and leg-to-leg fillet.
    """
    d: float
    b: float
    t: float
    s: float
    name: str = None

    # @todo - other properties to add: 
    # From AISC database: ro, H

    @classmethod
    def from_name(cls, name: str):
        key = name.upper()
        data = database.double_angle_database[key]
        return cls(data['d'], data['b'], data['t'], data['s'], name=key)

    @property
    def A(self):
        A = 2*(self.d+self.b-self.t)*self.t
        return A

    @property
    def y_bar(self):
        y_bar = (((self.b)*(self.t**2))-(self.t**3)+((self.d**2)*(self.t)))/(2*((self.b+self.d-self.t)*self.t))
        return y_bar

    @property
    def yp(self):
        if (self.b*self.t) > ((self.d-self.t)*self.t):
            yp = (0.5*(self.d+self.b-self.t)*self.t)/self.b
        else:
            yp = self.d - (0.5*(self.d+self.b-self.t)*self.t)/self.t
        return yp

    @property
    def Ix(self):
        Ix = 2*(((self.t/3)*((self.b*(self.t**2))+(self.d**3)-(self.t**3)))-((self.b+self.d-self.t)*self.t*(((((self.b)*(self.t**2))-(self.t**3)+((self.d**2)*(self.t)))/(2*((self.b+self.d-self.t)*self.t)))**2)))
        return Ix

    @property
    def Zx(self):
        if self.t <= ((self.d+self.b-self.t)*self.t/(2*self.b)):
            Zx = self.t*(((self.d-self.t)**2)-(self.b**2)+(2*self.b*self.d))/4
        else:
            Zx = (self.b*(self.t**2)/4)+((self.d*self.t*(self.d-self.t))/2)-(((self.t**2)*((self.d-self.t)**2))/(4*self.b))
        return 2*Zx

    @property
    def Sx(self):
        Sx = self.Ix/(self.d-self.y_bar)
        return Sx

    @property
    def rx(self):
        return sqrt(self.Ix/self.A)

    @property
    def Iy(self):
        Iy = (2*((((self.d-self.t)*(self.t**3))/12)+((self.d-self.t)*self.t*(((self.t/2)+(self.s/2))**2))))+(2*((((self.t)*(self.b**3))/12)+((self.t)*self.b*(((self.b/2)+(self.s/2))**2))))
        return Iy

    @property
    def Zy(self):
        Zy= 2*(((self.b-self.t)*self.t*((self.s/2)+self.b-((self.b-self.t)/2)))+(self.d*self.t*((self.s/2)+(self.t/2))))
        return Zy

    @property
    def Sy(self):
        Sy= (self.Iy)/((self.b+(self.s/2)))
        return Sy

    @property
    def ry(self):
        return sqrt(self.Iy/self.A)
