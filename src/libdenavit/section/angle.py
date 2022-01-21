import dataclasses
from math import sqrt,atan,radians,tan

from libdenavit.section import database


@dataclasses.dataclass
class Angle:
    """
    Parameters
    ----------
    d : float
        leg length along the x-axis
    b : float
        leg length along the y-axis
    t : float
        thickness
    name : str, optional
        name of the section

    Note that the definition of d and b are swapped from the DoubleAngle class.

    Calculations neglect the toe fillet and leg-to-leg fillet.
    """
    d: float
    b: float
    t: float
    name: str = None

    # @todo - other properties to add: 
    # From AISC database: Sz, ro, H, Iw

    @classmethod
    def from_name(cls, name: str):
        key = name.upper()
        data = database.angle_database[key]
        return cls(data['d'], data['b'], data['t'], name=key)

    @property
    def A(self):
        A = (self.d+self.b-self.t)*self.t
        return A

    @property
    def x_bar(self):
        x_bar =  (((self.b)*(self.t**2))-(self.t**3)+((self.d**2)*(self.t)))/(2*((self.d+self.b-self.t)*self.t))
        return x_bar

    @property
    def y_bar(self):
        y_bar = (((self.d)*(self.t**2))-(self.t**3)+((self.b**2)*(self.t)))/(2*((self.d+self.b-self.t)*self.t))
        return y_bar

    @property
    def xp(self):
        if self.b*self.t > 0.5*self.A:
            # PNA in vertical leg
            xp = 0.5*self.A/self.b
        else:
            xp = self.d - 0.5*self.A/self.t
        return xp

    @property
    def yp(self):
        if self.d*self.t > 0.5*self.A:
            # PNA in horizontal leg
            yp = 0.5*self.A/self.d
        else:
            yp = self.b - 0.5*self.A/self.t
        return yp

    @property
    def Ix(self):
        Ix = ((self.t/3)*((self.d*(self.t**2))+(self.b**3)-(self.t**3)))-((self.A)*(self.y_bar**2))
        return Ix

    @property
    def Zx(self):
        if self.t <= (self.A/(2*self.d)):
            Zx = self.t*(((self.b-self.t)**2)-(self.d**2)+(2*self.d*self.b))/4
        else:
            Zx = (self.d*(self.t**2)/4)+((self.b*self.t*(self.b-self.t))/2)-(((self.t**2)*((self.b-self.t)**2))/(4*self.d))
        return Zx

    @property
    def Sx(self):
        Sx = self.Ix/(self.b-self.y_bar)
        return Sx

    @property
    def rx(self):
        return sqrt(self.Ix/self.A)

    @property
    def Iy(self):
        Iy = ((self.t/3)*((self.b*(self.t**2))+(self.d**3)-(self.t**3)))-((self.A)*(self.x_bar**2))
        return Iy

    @property
    def Zy(self):
        if self.t <= (self.A/(2*self.b)):
            Zy= self.t*(((self.d-self.t)**2)-(self.b**2)+(2*self.b*self.d))/4
        else:
            Zy= (self.b*(self.t**2)/4)+((self.d*self.t*(self.d-self.t))/2)-(((self.t**2)*((self.d-self.t)**2))/(4*self.b))
        return Zy

    @property
    def Sy(self):
        Sy= self.Iy/(self.d-self.x_bar)
        return Sy

    @property
    def ry(self):
        return sqrt(self.Iy/self.A)

    @property
    def Ixy(self):
        Ixy = (((self.t**2)/4)*((self.d**2)+(self.b**2)-(self.t**2)))-(self.A*self.x_bar*self.y_bar)
        return Ixy

    @property
    def Iz(self):
        Iz = ((self.Ix+self.Iy)/2)-(sqrt((((self.Ix-self.Iy)/2)**2)+((abs(self.Ixy))**2)))
        return Iz

    @property
    def rz(self):
        return sqrt(self.Iz/self.A)

    @property
    def J(self):
        # CISC. (2002). Torsional Section Properties of Steel Shapes. Canadian Institute of Steel Construction, Ontario, Canada.
        d_d = self.d-(self.t/2)
        b_d = self.b-(self.t/2)
        J = (d_d+b_d)*(self.t**3)/3
        return J

    @property
    def Cw(self):
        Cw = ((self.t**3)/36)*((self.d-(self.t/2))**3+(self.b-(self.t/2))**3)
        return Cw

    @property
    def tan_alpha(self):
        if self.b == self.d:
            tan_alpha= tan(radians(45))
        else:
            tan_alpha= tan((atan((-2*self.Ixy)/(self.Ix-self.Iy)))/2)
        return tan_alpha


def compare_to_database():
    properties = ['A', 'x_bar', 'y_bar', 'xp', 'yp', 'Ix', 'Zx', 'Sx', 'rx', 'Iy', 'Zy', 'Sy', 'ry', 'Iz', 'rz', 'J', 'Cw', 'tan_alpha']

    for prop in properties:
        print('\n=== Checking %s ===' % prop)

        max_error_upper = 0.
        max_error_lower = 0.

        for key, iAngle in database.angle_database.items(): 
            
            # Get property from Python class
            s = Angle.from_name(key)
            X_calc = getattr(s, prop)
            
            # Get proproety from database
            if prop == 'x_bar':
                X_database = iAngle['x']
            elif prop == 'y_bar':
                X_database = iAngle['y']
            else:
                X_database = iAngle[prop]    
            
            # Compare
            percent_diff = 100*(X_calc-X_database)/X_database
            
            if abs(percent_diff) > 4:
                print(f'{key} --- {X_calc:.4f} / {X_database:.4f} --- {percent_diff:.4f}%')
            
            if percent_diff > max_error_upper:
                max_error_upper = percent_diff
            if percent_diff < max_error_lower:
                max_error_lower = percent_diff

        print('Error Summary:')
        print(f'Upper limit: {max_error_upper:.4f}%')
        print(f'Lower limit: {max_error_lower:.4f}%')

if __name__ == "__main__":
    compare_to_database()


