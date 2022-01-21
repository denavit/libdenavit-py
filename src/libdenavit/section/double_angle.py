import dataclasses
from math import pi,sqrt,atan,radians,tan
from libdenavit.section import Angle,database
from libdenavit.design import available_strength

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
    
    @property
    def rz_single(self):
        a = Angle(self.b,self.d,self.t)
        return a.rz
    
    @property
    def J(self):
        # CISC. (2002). Torsional Section Properties of Steel Shapes. Canadian Institute of Steel Construction, Ontario, Canada.
        d_d = self.d-(self.t/2)
        b_d = self.b-(self.t/2)
        J = 2*((d_d+b_d)*(self.t**3)/3)
        return J
    
    @property
    def yo(self):
        yo = self.y_bar- self.t/2
        return yo
    
    @property
    def ro(self):
        # Note xo = 0
        ro = sqrt(self.yo**2 + (self.Ix+self.Iy)/self.A)
        return ro
    
    @property
    def H(self):
        # Note xo = 0
        H = 1 - (self.yo**2/self.ro**2)
        return H
        
        
class DoubleAngleMember_SJI2020:
    def __init__(self,section,Fy,E,L,strength_type):
        self.section = section
        self.Fy = Fy
        self.E = E
        self.L = L
        self.strength_type = strength_type

    def Pnt(self):
        Pn = self.Fy*self.section.A
        return available_strength(Pn,self.strength_type,0.9,1.67)

    def Pnc(self):
        # Currently assumes that a sufficient number of fillers are provided
        # such that buckling of individual angle (rz) does not control
        # @todo - add buckling of individual angle     
        
        b_over_t = max(self.section.d,self.section.b)/self.section.t
        if b_over_t <= 0.45*sqrt(self.E/self.Fy):
            Q = 1.0
        elif b_over_t <= 0.91*sqrt(self.E/self.Fy):
            Q = 1.34-0.76*b_over_t*sqrt(self.Fy/self.E)
        else:
            Q = 0.53*self.E/(self.Fy*b_over_t**2)

        if self.L == 0.0:
            Fcr = Q*self.Fy
        else:
            KL_over_r = 1.0*self.L/self.section.rx
            Fe = pi**2*self.E/KL_over_r**2 
            if KL_over_r <= 4.71*sqrt(self.E/(Q*self.Fy)):
                Fcr = Q*(0.658**(Q*self.Fy/Fe))*self.Fy
            else:
                Fcr = 0.877*Fe
                
        Pn = Fcr*self.section.A
        return available_strength(Pn,self.strength_type,0.9,1.67)



def compare_to_database():
    properties = ['A', 'y_bar', 'yp', 'Ix', 'Zx', 'Sx', 'rx', 'Iy', 'Zy', 'Sy', 'ry', 'ro', 'H']

    for prop in properties:
        print('\n=== Checking %s ===' % prop)

        max_error_upper = 0.
        max_error_lower = 0.

        for key, iDoubleAngle in database.double_angle_database.items(): 
            
            # Get property from Python class
            s = DoubleAngle.from_name(key)
            X_calc = getattr(s, prop)

            # Get property from database
            if prop == 'y_bar':
                X_database = iDoubleAngle['y']
            else:
                X_database = iDoubleAngle[prop] 
            
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
