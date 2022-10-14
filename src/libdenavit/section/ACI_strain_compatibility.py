import numpy as np
from math import sin, cos, radians, pi


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
        for value in strain:
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
        for value in strain:
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
    max_compressive_strength = None

    def __init__(self, fiber_section, axes_origin="AsDefined"):
        self.fiber_section = fiber_section
        self.axes_origin = axes_origin

    def add_concrete_boundary(self, x, y, r):
        if isinstance(x, list):
            if len(x) == len(y) == len(r):
                for i in range(len(x)):
                    self._concrete_boundary_x = np.vstack((self._concrete_boundary_x, x[i]))
                    self._concrete_boundary_y = np.vstack((self._concrete_boundary_y, y[i]))
                    self._concrete_boundary_r = np.vstack((self._concrete_boundary_r, r[i]))
        else:
            self._concrete_boundary_x = np.vstack((self._concrete_boundary_x, x))
            self._concrete_boundary_y = np.vstack((self._concrete_boundary_y, y))
            self._concrete_boundary_r = np.vstack((self._concrete_boundary_r, r))

    def add_steel_boundary(self, x, y, r):
        self._steel_boundary_x = np.vstack((self._steel_boundary_x, x))
        self._steel_boundary_y = np.vstack((self._steel_boundary_y, y))
        self._steel_boundary_r = np.vstack((self._steel_boundary_r, r))

    def add_material(self, name, material_type, *args):
        if type(material_type).__name__ == "AciStrainCompatibilitySteelMaterial" or \
           type(material_type).__name__ == "AciStrainCompatibilityConcreteMaterial":
            mat = material_type
            
        elif material_type.lower() == 'steel':
            mat = AciStrainCompatibilitySteelMaterial(args[0], args[1])

        elif material_type.lower() == 'concrete':
            mat = AciStrainCompatibilityConcreteMaterial(args[0], args[1])

        else:
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

    def compute_point(self, xpt, ypt, angle, degrees=False):

        if degrees:
            angle = radians(angle)

        # Compute the strain of each fiber        
        a =  sin(angle)
        b = -cos(angle)
        c = -sin(angle) * xpt + cos(angle) * ypt
        y_ecf = self.extreme_concrete_compression_fiber(xpt, ypt, angle)
        if y_ecf < 0:
            strain = self.extreme_concrete_compression_strain / y_ecf * (a * self.x + b * self.y + c)
        else:
            strain = self.default_tensile_strain * np.ones(self.m.size)

        # Use the _materials to get stress for each fiber
        stress = np.zeros([len(self.m)])
        for i in self.uniq_mats:
            # Find fibers of the material
            ind = np.where(self.m == i)
            stress[ind] = self._materials[i].get_stress(strain[ind])

        # Sum to get P, Mx, and My
        P  = np.sum(stress * self.A)
        Mx = np.sum(stress * self.A * -self.y)
        My = np.sum(stress * self.A *  self.x)
        et = self.extreme_steel_tensile_strain(xpt, ypt, angle)

        return P, Mx, My, et

    def compute_point_uniform(self, strain):

        # Use the _materials to get stress for each fiber
        stress = np.zeros([len(self.m)])
        for i in self.uniq_mats:
            # Find fibers of the material
            ind = np.where(self.m == i)
            stress[ind] = self._materials[i].get_stress([strain])

        # Sum to get P, Mx, and My
        P  = np.sum(stress * self.A)
        Mx = np.sum(stress * self.A * -self.y)
        My = np.sum(stress * self.A *  self.x)
        et = strain

        return P, Mx, My, et
        
    def compute_section_interaction_2d(self,angle,num_points,degrees=False):

        if degrees:
            angle = radians(angle)

        # Initlize output
        P_list = []
        Mx_list = []
        My_list = []
        et_list = []

        # Select neutral axis locations
        ymin,ymax = self.fiber_section.get_bounds_at_angle(angle)
        d = ymax-ymin
        middle_points = np.linspace(ymin-0.51*d,ymax+0.51*d,num_points-2)
        # Uniform Compression
        P, Mx, My, et = self.compute_point_uniform(self.extreme_concrete_compression_strain)
        P_list.append(P)
        Mx_list.append(Mx)
        My_list.append(My)
        et_list.append(et)
                       
        # Middle Points
        for point in middle_points:
            xpt = -sin(angle)*point
            ypt =  cos(angle)*point
            P, Mx, My, et = self.compute_point(xpt, ypt, angle)
            P_list.append(P)
            Mx_list.append(Mx)
            My_list.append(My)
            et_list.append(et)        
        
        # Uniform Tension
        P, Mx, My, et = self.compute_point_uniform(self.default_tensile_strain)
        P_list.append(P)
        Mx_list.append(Mx)
        My_list.append(My)
        et_list.append(et)
        
        # Middle Points
        for point in middle_points:
            xpt = -sin(angle)*point
            ypt =  cos(angle)*point
            P, Mx, My, et = self.compute_point(xpt, ypt, angle+pi)
            P_list.append(P)
            Mx_list.append(Mx)
            My_list.append(My)
            et_list.append(et)  

        # Uniform Compression
        P, Mx, My, et = self.compute_point_uniform(self.extreme_concrete_compression_strain)
        P_list.append(P)
        Mx_list.append(Mx)
        My_list.append(My)
        et_list.append(et)
        
        # Convert lists to numpy arrays
        P_list = np.array(P_list)
        Mx_list = np.array(Mx_list)
        My_list = np.array(My_list)
        et_list = np.array(et_list)
        
        # Apply Maximum Compressive Strength Cap
        if self.max_compressive_strength is not None:
            if self.max_compressive_strength > 0:
                raise ValueError(f"max_compressive_strength set to a positive (i.e., tensile) value: {self.max_compressive_strength}")
            for i, v in enumerate(P_list):
                if v < self.max_compressive_strength:
                    P_list[i] = self.max_compressive_strength

        return P_list, Mx_list, My_list, et_list