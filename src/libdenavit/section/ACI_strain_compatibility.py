import numpy as np
from math import sin, cos


class AciStrainCompatibilitySteelMaterial:

    def __init__(self, Fy, Es):
        self.Fy = Fy  # float
        self.Es = Es  # float

    @property
    def ey(self):
        # Returns yield strain
        return self.Fy / self.Es

    def get_stress(self, strain):
        # Input is an array of strains
        # Output is an array of stresses
        # Follow an elastic perfectly plastic material model

        stress = []
        for index, value in enumerate(strain):
            if value <= -self.ey:
                stress.append(-self.Fy)
            elif value <= self.ey:
                stress.append(value * self.Es)
            else:
                stress.append(self.Fy)

        return stress


class AciStrainCompatibilityConcreteMaterial:

    def __init__(self, fc, units):
        self.fc = fc  # float
        self.units = units  # string

    extreme_concrete_compression_strain = -0.003

    @property
    def beta1(self):
        # Returns beta1 per ACI 318-19 Table 22.2.2.4.3 with the exception that 
        # beta1 is taken as 0.85 for fc < 2500 psi instead of being undefined. 
        
        if self.units.lower() == "us":
            # fc in units of ksi
            if self.fc <= 4:
                beta1 = 0.85
            elif self.fc <= 8:
                beta1 = 0.85 - 0.05 * (self.fc - 4)
            else:
                beta1 = 0.65

        elif self.units.lower() == "si":
            # fc in units of MPa
            if self.fc <= 28:
                beta1 = 0.85
            elif self.fc <= 55:
                beta1 = 0.85 - 0.05 * (self.fc - 28) / 7
            else:
                beta1 = 0.65

        else:
            raise ValueError(f"The given units ({self.units}) is not supported yet")

        return beta1

    def get_stress(self, strain):
        # Input is an array of strains
        # Output is an array of equivalent stresses based on concrete stress block
        
        # Strain at which the concrete stress block initiates 
        ecr = self.extreme_concrete_compression_strain * (1 - self.beta1)

        stress = []
        for index, value in enumerate(strain):
            if value <= ecr:
                stress.append(-0.85 * self.fc)
            else:
                stress.append(0)

        return stress


class AciStrainCompatibility:
    # Several functions in this class take xpt, ypt, angle as input xpt, ypt, angle are floats xpt and ypt are
    # the coordinates of any point on the neutral axis angle is the angle in radians between the positive x axis and
    # the neutral axis. in the direction of the angle, compression is to the left and tension is to the right For
    # example angle = 0 indicates bending about the x-axis (or a line parallel to the x-axis) with compression above
    # the neutral axis and tension below the neutral axis

    _concrete_boundary_x = np.zeros((0, 1))
    _concrete_boundary_y = np.zeros((0, 1))
    _concrete_boundary_r = np.zeros((0, 1))

    _steel_boundary_x = np.zeros((0, 1))
    _steel_boundary_y = np.zeros((0, 1))
    _steel_boundary_r = np.zeros((0, 1))

    _materials = dict()

    extreme_concrete_compression_strain = -0.003
    default_tensile_strain = 0.005

    def __init__(self, fiber_section, axes_origin="AsDefined"):
        self.fiber_section = fiber_section
        self.axes_origin = axes_origin

    def add_concrete_boundary(self, x, y, r):
        self._concrete_boundary_x = np.vstack((self._concrete_boundary_x, x))
        self._concrete_boundary_y = np.vstack((self._concrete_boundary_y, y))
        self._concrete_boundary_r = np.vstack((self._concrete_boundary_r, r))

    def add_steel_boundary(self, x, y, r):
        self._steel_boundary_x = np.vstack((self._steel_boundary_x, x))
        self._steel_boundary_y = np.vstack((self._steel_boundary_y, y))
        self._steel_boundary_r = np.vstack((self._steel_boundary_r, r))

    def add_material(self, name, material_type, *args):
        if type(material_type) == "ACI_strain_compatibility.AciStrainCompatibilitySteelMaterial" or \
                type(material_type) == "ACI_strain_compatibility.AciStrainCompatibilityConcreteMaterial":
            mat = material_type
            
        elif material_type.lower() == 'steel':
            mat = AciStrainCompatibilitySteelMaterial(args[0], args[1])

        elif material_type.lower() == 'concrete':
            mat = AciStrainCompatibilityConcreteMaterial(args[0], args[1])

        else:
            # Error unknown material type
            raise ValueError(f"Unknown material type: {material_type}")

        self._materials[name] = mat

    def extreme_concrete_compression_fiber(self, xpt, ypt, angle):
        # Returns the distance from the neutral axis to the extreme concrete compression fiber
        # Despite the name, it is not based on the fiber locations, it is based on the true 
        # edge of material, which is defined by the concrete boundary points

        a =  sin(angle)
        b = -cos(angle)
        c = -sin(angle) * xpt + cos(angle) * ypt

        y = a * self._concrete_boundary_x + b * self._concrete_boundary_y + c
        y = np.amin(y - self._concrete_boundary_r)
        return y

    def extreme_steel_tension_fiber(self, xpt, ypt, angle):
        # Returns the distance from the neutral axis to the extreme steel tension fiber
        # Despite the name, it is not based on the fiber locations, it is based on the true 
        # edge of material, which is defined by the steel boundary points   

        a =  sin(angle)
        b = -cos(angle)
        c = -sin(angle) * xpt + cos(angle) * ypt

        y = a * self._steel_boundary_x + b * self._steel_boundary_y + c
        y = np.amax(y + self._steel_boundary_r)
        return y

    def extreme_steel_tensile_strain(self, xpt, ypt, angle):
        # Returns the strain in the extreme tension fiber when the extreme concrete compression strain equals the
        # defined value If the neutral axis is defined such that no concrete is in compression, then it returns
        # self.default_tensile_strain

        yc = self.extreme_concrete_compression_fiber(xpt, ypt, angle)
        yt = self.extreme_steel_tension_fiber(xpt, ypt, angle)
        ec = self.extreme_concrete_compression_strain
        if yc < 0:
            et = yt * (ec / yc)
        else:
            et = self.default_tensile_strain

        return et

    def build_data(self):
        self.uniq_mats = self.fiber_section.unique_mat_ids()
        (self.A, self.x, self.y, self.m) = self.fiber_section.get_fiber_data()
        return

    def compute_point(self, xpt, ypt, angle):

        # Get fiber data
        # Compute perpendicular distance of each fiber from the neutral axis
        # Compute the strain of each fiber
        # Use the _materials to get stress for each fiber
        # Sum to get P, Mx, and My
        # P, Mx, and My are floats
        a =  sin(angle)
        b = -cos(angle)
        c = -sin(angle) * xpt + cos(angle) * ypt
        y_ecf = self.extreme_concrete_compression_fiber(xpt, ypt, angle)
        if y_ecf < 0:
            strain = self.extreme_concrete_compression_strain / y_ecf * (a * self.x + b * self.y + c)
        else:
            strain = self.default_tensile_strain * np.ones(self.m.size)

        stress = np.zeros([len(self.m)])
        for i in self.uniq_mats:
            # Find fibers of the material
            ind = np.where(self.m == i)
            stress[ind] = self._materials[i].get_stress(strain[ind])

        P  = sum(stress * self.A)
        Mx = sum(stress * self.A * - self.y)
        My = sum(stress * self.A * self.x)
        et = self.extreme_steel_tensile_strain(xpt, ypt, angle)

        return P, Mx, My, et


def run_example():
    import matplotlib.pyplot as plt
    from libdenavit.section import FiberQuadPatch, FiberSingle, FiberSection

    # Parameter
    B = 20 # Cross-sectional width
    H = 40 # Cross-sectional height
    fc = 4
    fy = 60
    Es = 29000
    
    # Build Fiber Section Object
    c = FiberSection(50, 50)
    b1 = FiberQuadPatch(-10, -20, -10, 20, 10, 20, 10, -20, 1)  # xI, yI, ..., m
    a1 = FiberSingle(12, -4, -14, 2, 1)  # A, x, y, m, m_neg
    a2 = FiberSingle(12, -4,  14, 2, 1)
    a3 = FiberSingle(12,  4, -14, 2, 1)
    a4 = FiberSingle(12,  4,  14, 2, 1)
    c.add_fibers(b1, a1, a2, a3, a4)

    # Output Information From Fiber Section
    c.print_section_properties()
    c.plot_fibers()

    # Build ACI Strain Compatibility Object
    test = AciStrainCompatibility(c)
    test.add_concrete_boundary( 0.5*B, 0.5*H, 0)
    test.add_concrete_boundary( 0.5*B,-0.5*H, 0)
    test.add_concrete_boundary(-0.5*B,-0.5*H, 0)
    test.add_concrete_boundary(-0.5*B, 0.5*H, 0)
    test.add_steel_boundary(4, -14, 0)
    test.add_steel_boundary(4, 14, 0)
    test.add_steel_boundary(-4, 14, 0)
    test.add_steel_boundary(-4, -14, 0)
    test.add_material(1, "concrete", fc, "US")
    test.add_material(2, "steel", fy, Es)
    test.build_data()
    
    # Select Neutral Axis Locations
    a = np.arange(-70, -30, 0.5)
    b = np.arange(-30,  20, 0.1)
    c = np.arange( 20,  80, 0.1)
    d = np.concatenate((a, b, c))
    
    # Perform Interaction Calculations
    plot_P = []
    plot_Mx = []  
    for i in d:
        P, Mx, My, et = test.compute_point(0, i, 0)
        plot_P.append(-P)
        plot_Mx.append(Mx)

    # Plot Interaction Diagram
    plt.plot(plot_Mx, plot_P, 'bo-')
    plt.show()


if __name__ == "__main__":
    run_example()


