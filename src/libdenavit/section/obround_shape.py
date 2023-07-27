import numpy as np
from libdenavit.section import GeometricShape, FiberCirclePatch, FiberQuadPatch
import matplotlib.pyplot as plt


class Obround(GeometricShape):
    """
    Parameters
    ----------
    D : float
        diameter of semicircles
    a : float
        top and bottom length of interior rectangle
    """
    
    section_type = "Obround"
    
    def __init__(self, D, a):
        if D <= 0:
            raise ValueError("D must be greater than 0")
        if a <= 0:
            raise ValueError("a must be greater than 0")
        self.D = D
        self.a = a
    
    @property
    def is_section_valid(self):
        tf1 = self.D > 0  # diameter should be positive
        tf2 = self.a > 0  # length of flat should be positive
        tf = all([tf1,tf2])
        return tf
    
    def depth(self, axis):
        if axis.lower() == 'x':
            return self.D
        elif axis.lower() == 'y':
            return self.D + self.a
        else:
            raise ValueError(f'Unknown axis: {axis}')

    @property
    def perimeter(self):
        p = 2 * (np.pi * self.D/2 + self.a)
        return p
    @property
    def A(self):
        a = (np.pi / 4) * self.D ** 2 + self.D*self.a
        return a
    
    @property
    def Ix(self):
        i = (np.pi / 64) * self.D ** 4 + (1/12) * self.a * self.D ** 3
        return i
        
    @property
    def Iy(self):
        # each semicircle
        i1 = self.D**4*(np.pi/128 - 1/(18*np.pi)) + (np.pi/8)*self.D**2 * (self.a/2 + (2*self.D)/(3*np.pi))**2
        # central rectangle
        i2 = (1/12)*self.D*self.a**3
        # whole shape
        i = 2*i1 + i2
        return i    
    
    @property
    def boundary_points(self):
        x = [-self.a/2,self.a/2]
        y = [0,0]
        r = [self.D/2,self.D/2]
        return x, y, r
    
    def plot_section(self, **kwargs):
        r = self.D / 2
        theta = np.linspace(-np.pi/2, np.pi/2, 40)
        x1 =  self.a/2 + r * np.cos(theta)
        y1 = r * np.sin(theta)
        x2 = -self.a/2 + r * np.cos(theta+np.pi)
        y2 = r * np.sin(theta+np.pi)        
        x = np.append(x1,x2)
        y = np.append(y1,y2)
        plt.fill(x, y, **kwargs)
    
    def add_to_fiber_section(self, fiber_section, mat_id):
        f1 = FiberCirclePatch(-self.a/2, 0, 0, self.D/2, mat_id, a1= np.pi/2, a2=3*np.pi/2)
        f2 = FiberCirclePatch( self.a/2, 0, 0, self.D/2, mat_id, a1=-np.pi/2, a2= np.pi/2)
        f3 = FiberQuadPatch(-self.a/2, -self.D/2, -self.a/2, self.D/2, self.a/2, self.D/2, self.a/2, -self.D/2, mat_id)
        fiber_section.add_fibers(f1,f2,f3)


if __name__ == "__main__":
    shape = Obround(24,12)
    shape.plot_section(facecolor='lightsalmon', edgecolor='orangered', linewidth=3)
    plt.show()
