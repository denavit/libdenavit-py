import numpy as np
from libdenavit.section import GeometricShape, FiberCirclePatch
import matplotlib.pyplot as plt


class Circle(GeometricShape):
    """
    Parameters
    ----------
    D : float
        diameter of circle
    """
    
    section_type = "Circle"
    
    def __init__(self, D):
        if D <= 0:
            raise ValueError("Diameter must be greater than 0")
        self.diameter = D
    
    @property
    def is_section_valid(self):
        tf1 = self.diameter > 0  # diameter should be positive
        tf = all([tf1])
        return tf
    
    def depth(self, axis):
        return self.diameter

    @property
    def perimeter(self):
        p = np.pi * self.diameter
        return p

    @property
    def A(self):
        a = (np.pi / 4) * self.diameter ** 2
        return a
    
    @property
    def Ix(self):
        i = (np.pi / 64) * self.diameter ** 4
        return i
        
    @property
    def Iy(self):
        i = (np.pi / 64) * self.diameter ** 4
        return i    
    
    @property
    def boundary_points(self):
        x = 0
        y = 0
        r = self.diameter / 2
        return x, y, r
    
    def plot_section(self, **kwargs):
        theta = np.linspace(0, 2 * np.pi, 80)
        r = self.diameter / 2
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        plt.fill(x, y, **kwargs)
    
    def add_to_fiber_section(self, fiber_section, mat_id):
        f = FiberCirclePatch(0, 0, 0, self.diameter/2, mat_id)
        fiber_section.add_fibers(f)


if __name__ == "__main__":
    circle = Circle(2)
    circle.plot_section(facecolor='lightsalmon', edgecolor='orangered', linewidth=3)
    plt.show()
