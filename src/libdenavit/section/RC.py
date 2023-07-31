from math import sqrt, pi, ceil, exp, sin
import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops
from libdenavit.OpenSees import circ_patch_2d, obround_patch_2d, obround_patch_2d_confined
from libdenavit.section import AciStrainCompatibility, FiberSection, ACI_phi


class RC:
    treat_reinforcement_as_point = True
    _reinforcement = None
    _Ec = None
    _Es = None
    _eps_c = None
    _Abt = None

    def __init__(self, conc_cross_section, reinforcement, fc, fy, units, dbt=None, s=None, fyt=None, lat_config="A",
                 transverse_reinf_type='ties'):
        self.conc_cross_section = conc_cross_section
        self.reinforcement = reinforcement
        self.fc = fc
        self.fy = fy
        self.units = units
        self.dbt = dbt
        self.s = s
        self.fyt = fyt
        self.lat_config = lat_config
        self.transverse_reinf_type = transverse_reinf_type

    @property
    def Ec(self):
        if self._Ec is not None:
            return self._Ec

        if self.units.lower() == 'us':
            return 57 * sqrt(self.fc * 1000)

        if self.units.lower() == 'si':
            return 4700 * sqrt(self.fc)

        raise ValueError(f'Ec is not set and default is not implemented for {self.units = }')

    @Ec.setter
    def Ec(self, x):
        self._Ec = x

    @property
    def Abt(self):
        if self._Abt is not None:
            return self._Abt
        else:
            return pi * self.dbt ** 2 / 4
    
    @Abt.setter
    def Abt(self, x):
        self._Abt = x
    
    @property
    def Es(self):
        if self._Es is not None:
            return self._Es

        if self.units.lower() == 'us':
            return 29000

        if self.units.lower() == 'si':
            return 2e5

        raise ValueError(f'Es is not set and default is not implemented for {self.units = }')

    @Es.setter
    def Es(self, x):
        self._Es = x
    
    @property
    def eps_c(self):
        if self._eps_c is not None:
            return self._eps_c
        if self.units.lower() == 'us':
            return (self.fc * 1000) ** (1 / 4) / 4000
        if self.units.lower() == 'si':
            return self.fc ** (1 / 4) / 28

        raise ValueError(f'eps_c is not set and default is not impleted for {self.units = }')
    
    @eps_c.setter
    def eps_c(self, x):
        self._eps_c = x
        
    @property
    def reinforcement(self):
        return self._reinforcement

    @reinforcement.setter
    def reinforcement(self, x):
        if type(x) == list:
            self._reinforcement = x
        else:
            self._reinforcement = [x]

    def depth(self, axis):
        d = self.conc_cross_section.depth(axis)
        return d

    @property
    def Ag(self):
        return self.conc_cross_section.A

    @property
    def Ac(self):
        return self.Ag - self.Asr

    @property
    def Asr(self):
        a = 0
        for i in self.reinforcement:
            a += i.num_bars * i.Ab
        return a
    def maximum_concrete_compression_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        if type(self.conc_cross_section).__name__ == 'Rectangle':
            extreme_strain = axial_strain - self.conc_cross_section.H/2 * abs(curvatureX) \
                             - self.conc_cross_section.B/2 * abs(curvatureY)
            return extreme_strain

        if type(self.conc_cross_section).__name__ == 'Circle':
            extreme_strain = axial_strain - self.conc_cross_section.diameter/2 * sqrt(curvatureX**2 + curvatureY**2)
            return extreme_strain

        if  type(self.conc_cross_section).__name__ == 'Obround':
            extreme_strain = float("inf")
            for i in range(100):
                x_coord = self.conc_cross_section.a/2 + self.conc_cross_section.D/2 * sin(pi/2 * i/100)
                y_coord = self.conc_cross_section.D/2 - self.conc_cross_section.D/2 * sin(pi/2 * i/100)
                strain = axial_strain - y_coord * abs(curvatureX) - x_coord * abs(curvatureY)
                if strain < extreme_strain:
                    extreme_strain = strain
            return extreme_strain

        raise ValueError(f'No maximum concrete compression strain implemented for {type(self.conc_cross_section).__name__ = }')

    def maximum_tensile_steel_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        max_strain = float('-inf')
        for i in self.reinforcement[0].coordinates:
            strain = axial_strain - (i[1] - self.reinforcement[0].yc) * curvatureX \
                     - (i[0] - self.reinforcement[0].xc) * curvatureY
            if strain > max_strain:
                max_strain = strain
        return max_strain

    def Ig(self, axis):
        return self.conc_cross_section.I(axis)

    def Ic(self, axis):
        return self.Ig(axis) - self.Isr(axis)

    def Isr(self, axis):
        i = 0
        for j in self.reinforcement:
            i += j.I(axis)
        return i

    @property
    def p0(self):
        # Nominal axial strength of section
        p0 = 0.85 * self.fc * (self.Ag - self.Asr) + self.fy * self.Asr
        return p0

    @property
    def p0g(self):
        # Nominal axial strength of gross concrete section
        p0g = 0.85 * self.fc * self.Ag
        return p0g

    @property
    def pnco(self):
        # See Section 22.4.2 of ACI318-19
        if self.transverse_reinf_type.lower() == 'ties':
            pnco = 0.80 * self.p0
        elif self.transverse_reinf_type.lower() in ['spiral', 'spirals']:
            pnco = 0.85 * self.p0
        else:
            raise ValueError("Unknown transverse_reinf_type")
        return pnco

    def EIgross(self, axis):
        return self.Ec * self.Ic(axis) + self.Es * self.Isr(axis)
    
    def EIeff(self, axis, EI_type, betadns=0.0, P=None, M=None):
        if EI_type == "ACI-a":
            # ACI 318-19, Section 6.6.4.4.4a
            return 0.4 * self.Ec * self.Ig(axis) / (1 + betadns)

        if EI_type == "ACI-b":
            # ACI 318-19, Section 6.6.4.4.4b
            return (0.2 * self.Ec * self.Ig(axis) + self.Es * self.Isr(axis)) / (1 + betadns)

        if EI_type == "ACI-c":
            # ACI 318-19, Section 6.6.4.4.4c
            if P is None or M is None:
                raise ValueError("P and M must be defined for EI_type = 'ACI-c'")
            EI = []
            max_EI = 0.875 * self.Ec * self.Ig(axis)
            min_EI = 0.35 * self.Ec * self.Ig(axis)
            if isinstance(P, (list, np.ndarray)):
                if type(P) != type(M) or len(P) != len(M):
                    raise ValueError("P and M must be of the same type and size")
                for P, M in zip(P, M):
                        EI_ACI = (0.8+25*self.Asr/self.Ag)*(1-M/P/self.depth(axis)-0.5*P/self.p0)*self.Ec*self.Ig(axis)
                        if EI_ACI > max_EI:
                            EI.append(max_EI)
                        elif EI_ACI < min_EI:
                            EI.append(min_EI)
                        else:
                            EI.append(EI_ACI)
                return EI

            elif isinstance(P, float):
                EI_ACI = (0.8 + 25 * self.Asr / self.Ag) * (
                            1 - M / P / self.depth(axis) - 0.5 * P / self.p0) * self.Ec * self.Ig(axis)
                if EI_ACI > max_EI:
                    EI = max_EI
                elif EI_ACI < min_EI:
                    EI = min_EI
                else:
                    EI = EI_ACI
                return EI

        if EI_type == "JF-a":
            # Jenkins and Frosch, 2011 - Eq.10.1
            EI = []
            if P is None or M is None:
                raise ValueError("P and M must be defined for EI_type = 'JF-a'")
            elif isinstance(P, (list, np.ndarray)) and isinstance(M, (list, np.ndarray)):
                if len(P) != len(M):
                    raise ValueError("P and M must be of the same size")
                for P, M in zip(P, M):
                    ecc_ratio = abs(M) / abs(P) / self.depth(axis)
                    if ecc_ratio <= 0.1:
                        EI_jenkins = (1.05 - 0.6 * abs(P) / self.p0) * (1 + 3 * (self.Asr / self.Ag - 0.01)) *\
                                     self.Ec * self.Ig(axis) / (1 + betadns)
                        min_EI = 0.3 * self.Ec * self.Ig(axis) / (1 + betadns)
                        EI.append(max(EI_jenkins, min_EI))
                    else:
                        EI_jenkins = (1.05 - 0.6 * abs(P) / self.p0) * (1 + 3 * (self.Asr / self.Ag - 0.01)) * (
                                1.2 - 2 * ecc_ratio) * self.Ec * self.Ig(axis) / (1 + betadns)
                        min_EI = 0.3 * self.Ec * self.Ig(axis) / (1 + betadns)
                        EI.append(max(EI_jenkins, min_EI))
                return EI

            elif isinstance(P, (float, np.float64, int)) and isinstance(M, (float, np.float64, int)):
                try:
                    ecc_ratio = abs(M) / abs(P) / self.depth(axis)
                except ZeroDivisionError:
                    return np.nan
                if ecc_ratio <= 0.1:
                    EI_jenkins = (1.05 - 0.6 * abs(P) / self.p0) * (1 + 3 * (self.Asr / self.Ag - 0.01)) * \
                                 self.Ec * self.Ig(axis) / (1 + betadns)
                    min_EI = 0.3 * self.Ec * self.Ig(axis) / (1 + betadns)
                    EI = max(EI_jenkins, min_EI)
                else:
                    EI_jenkins = (1.05 - 0.6 * abs(P) / self.p0) * (1 + 3 * (self.Asr / self.Ag - 0.01)) * (
                            1.2 - 2 * ecc_ratio) * self.Ec * self.Ig(axis) / (1 + betadns)
                    min_EI = 0.3 * self.Ec * self.Ig(axis) / (1 + betadns)
                    EI = max(EI_jenkins, min_EI)
                return EI

            else:
                raise ValueError("P and M types or sizes are not as expected")

        if EI_type == "JF-b":
            # Jenkins and Frosch, 2011 - Eq.10.2
            EI = []
            if P is None or M is None:
                raise ValueError("P and M must be defined for EI_type = 'JF-b'")

            elif isinstance(P, (list, np.ndarray)) and isinstance(M, (list, np.ndarray)):
                if len(P) != len(M):
                    raise ValueError("P and M must be of the same size")
                for P, M in zip(P, M):
                    ecc_ratio = M / P / self.depth(axis)
                    if ecc_ratio <= 0.1:
                        EI_jenkins = (1.0 - 0.5 * P / self.p0) * self.Ec * self.Ig(axis) / (1 + betadns)
                        min_EI = 0.4 * self.Ec * self.Ig(axis) / (1 + betadns)
                        EI.append(max(EI_jenkins, min_EI))
                    else:
                        EI_jenkins = (1.0 - 0.5 * P / self.p0) * self.Ec * self.Ig(axis) * (
                            1.2 - 2 * ecc_ratio) / (1 + betadns)
                        min_EI = 0.4 * self.Ec * self.Ig(axis) / (1 + betadns)
                        EI.append(max(EI_jenkins, min_EI))
                return EI

            elif isinstance(P, (float, np.float64, int)) and isinstance(M, (float, np.float64, int)):
                try:
                    ecc_ratio = abs(M) / abs(P) / self.depth(axis)
                except ZeroDivisionError:
                    return np.nan
                if ecc_ratio <= 0.1:
                    EI_jenkins = (1.0 - 0.5 * P / self.p0) * self.Ec * self.Ig(axis) / (1 + betadns)
                    min_EI = 0.4 * self.Ec * self.Ig(axis) / (1 + betadns)
                    EI = max(EI_jenkins, min_EI)
                else:
                    EI_jenkins = (1.0 - 0.5 * P / self.p0) * self.Ec * self.Ig(axis) * (
                            1.2 - 2 * ecc_ratio) / (1 + betadns)
                    min_EI = 0.4 * self.Ec * self.Ig(axis) / (1 + betadns)
                    EI = max(EI_jenkins, min_EI)
                return EI

            else:
                raise ValueError("P and M types or sizes are not as expected")

        if EI_type == "gross":
            return self.EIgross(axis)

        raise ValueError(f'Unknown EI_type {EI_type}')

    def phi(self, et):
        f = ACI_phi(self.transverse_reinf_type, et, self.fy / self.Es)
        return f

    def plot_section(self, show=True, **kwargs):
        plt.figure()
        self.conc_cross_section.plot_section(edgecolor='k',facecolor=[0.9,0.9,0.9],**kwargs)
        for i in range(len(self.reinforcement)):
            self.reinforcement[i].plot_section(color='k',**kwargs)
        plt.box(False)
        plt.axis('scaled')
        if show:
            plt.show()

    def aci_strain_compatibility_object(self):
        id_conc = 1
        id_reinf = 2
        fs = self.fiber_section_object(id_conc, id_reinf)
        scACI = AciStrainCompatibility(fs)
        # add concrete boundaries
        x, y, r = self.conc_cross_section.boundary_points
        scACI.add_concrete_boundary(x, y, r)
        # add steel boundary
        for i in self.reinforcement:
            x, y = i.coordinates

            for j in range(len(x)):
                if self.treat_reinforcement_as_point:
                    scACI.add_steel_boundary(x[j], y[j], 0)
                else:
                    raise ValueError("Not implemented yet")

        scACI.max_compressive_strength = -self.pnco
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.fy, self.Es)
        return scACI

    def fiber_section_object(self, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        self.conc_cross_section.add_to_fiber_section(fs, id_conc)
        for i in self.reinforcement:
            i.add_to_fiber_section(fs, id_reinf, id_conc)
        return fs

    def section_interaction_2d(self,axis,num_points,factored=False, only_compressive=False):
        scACI = self.aci_strain_compatibility_object()
        scACI.build_data()
        
        if axis == 'x':
            P, M, _, et = scACI.compute_section_interaction_2d(0,num_points,degrees=True)
        elif axis == 'y':
            P, _, M, et = scACI.compute_section_interaction_2d(90,num_points,degrees=True)
        else:
            raise ValueError(f'Unknown axis ({axis})')
            
        if factored:
            ϕ = self.phi(et)
            P = ϕ*P
            M = ϕ*M

        if only_compressive:
            from libdenavit import InteractionDiagram2d
            P = -1 * P
            P_M_id2d = InteractionDiagram2d(M, P, is_closed=True)
            P_et_id2d = InteractionDiagram2d(et, P, is_closed=True)
            M_et_id2d = InteractionDiagram2d(M, et, is_closed=True)
            P = P[:-1]
            M = M[:-1]
            et = et[:-1]
            from libdenavit import find_limit_point_in_list
            ind, x = find_limit_point_in_list(P, 0)
            P = np.insert(P, ind+1, 0)
            M = np.insert(M, ind+1, P_M_id2d.find_x_given_y(0, 'pos'))
            et = np.insert(et, ind+1, P_et_id2d.find_x_given_y(0, 'pos')) # @todo- chcek this

            P = np.insert(P, 0, P_M_id2d.find_y_given_x(0, 'pos'))
            M = np.insert(M, 0, 0)
            et = np.insert(et, 0, M_et_id2d.find_y_given_x(0, 'pos'))# @todo- chcek this

            # Delete extra points
            ind_P_negative = np.where(P < 0)
            P = np.delete(P, ind_P_negative)
            M = np.delete(M, ind_P_negative)
            et = np.delete(et, ind_P_negative)
            ind_M_negative = np.where(M < 0)
            P = np.delete(P, ind_M_negative)
            M = np.delete(M, ind_M_negative)
            et = np.delete(et, ind_M_negative)
        return P, M, et

    def build_ops_fiber_section(self, section_id, start_material_id, steel_mat_type, conc_mat_type, nfy, nfx, GJ=1.0e6, axis=None):
        """ Builds the fiber section object

        Parameters
        ----------
        section_id : int
            The id of the section
        start_material_id : int
            The id of the first uniaxial material to be defined (others will be defined sequentially)
        steel_mat_type : str
            The type of the steel material
        conc_mat_type : str
            The type of the concrete material
        nfy : int
            The minimum number of fibers in the y direction
        nfx : int
            The minimum number of fibers in the x direction
        GJ : float (default 1.0e6)
            The torsional rigidity of the cross section
        axis : str (default None)
            Optional argument for defining a fiber section for 2-dimensional analysis
              - If "None", then a 3-dimensional fiber section will be defined
              - If "x", then a 2-dimensional fiber section will be defined based on 
                bending about the x-axis (note that the value nfx will be ignored) 
              - If "y", then a 2-dimensional fiber section will be defined based on 
                bending about the y-axis (note that the value nfy will be ignored)
        """

        # Two or three uniaxial materials are defined in this function
        steel_material_id = start_material_id
        concrete_material_id = start_material_id+1        # Used if no confinement
        cover_concrete_material_id = start_material_id+1  # Used if confinement
        core_concrete_material_id  = start_material_id+2  # Used if confinement
        
        # Define section        
        ops.section('Fiber', section_id, '-GJ', GJ)

        # region Define Steel Material
        if steel_mat_type == "ElasticPP":
            ops.uniaxialMaterial("ElasticPP", steel_material_id, self.Es, self.fy / self.Es)

        elif steel_mat_type == "Hardening":
            ops.uniaxialMaterial("Hardening", steel_material_id, self.Es, self.fy, 0.001*self.Es, 0)

        elif steel_mat_type == "ReinforcingSteel":
            ops.uniaxialMaterial("ReinforcingSteel", steel_material_id, self.fy, self.fy * 1.5, self.Es, self.Es / 2,
                                 0.002, 0.008)

        elif steel_mat_type == "Elastic":
            ops.uniaxialMaterial("Elastic", steel_material_id, self.Es)
        else:
            raise ValueError(f"Steel material {steel_mat_type} not supported")
        # endregion

        if type(self.conc_cross_section).__name__ == 'Rectangle':
            
            # region Define Concrete Material
            if conc_mat_type == "Concrete04":
                # Defined based on Mander, J. B., Priestley, M. J. N., and Park, R. (1988). 
                # “Theoretical Stress-Strain Model for Confined Concrete.” Journal of Structural 
                # Engineering, ASCE, 114(8), 1804―1826.
            
                if type(self.reinforcement[0]).__name__ != 'ReinfRect':
                    raise ValueError(f"Reinforcement type {type(self.reinforcement).__name__} not supported for this section type")
                if self.reinforcement[0].xc != 0 or self.reinforcement[0].yc != 0:
                    raise ValueError(f"Reinforcing pattern must be centered")
                if self.dbt is None:
                   raise ValueError("dbt must be defined")
                if self.s is None:
                    raise ValueError("s must be defined")
                
                # dc and bc = core dimensions to centerlines of perimeter hoop in x and y directions
                bc = self.reinforcement[0].Bx + self.reinforcement[0].db/2 + self.dbt/2
                dc = self.reinforcement[0].By + self.reinforcement[0].db/2 + self.dbt/2
                                
                # Ac = area of core of section enclosed by the center lines of the permiter hoops
                Ac = bc*dc
                
                # ρcc = ratio of area of longitudinal reinforcement to area of core of section
                ρcc = self.reinforcement[0].Ab * self.reinforcement[0].num_bars / Ac
                
                # wx and wy = clear distance between bars in x and y directions
                wx = self.reinforcement[0].Bx / (self.reinforcement[0].nbx - 1) - self.reinforcement[0].db
                wy = self.reinforcement[0].By / (self.reinforcement[0].nby - 1) - self.reinforcement[0].db
                
                # sp = clear vertical spacing between spiral or hoop bars
                sp = self.s - self.dbt
                
                # ke = confinement effectiveness coefficient
                sum_w = (2 * (self.reinforcement[0].nbx - 1) * wx ** 2 +
                         2 * (self.reinforcement[0].nby - 1) * wy ** 2) / (6 * bc * dc)
                ke = (1 - sum_w) * (1 - sp / (2 * bc)) * (1 - sp / (2 * dc)) / (1 - ρcc)
                
                # Asx and Asy = total area of transverse bars running in the x and y directions
                if self.lat_config == 'A':
                    Asx = 2 * pi * self.dbt ** 2 / 4
                    Asy = 2 * pi * self.dbt ** 2 / 4
                elif self.lat_config == 'B':
                    Asx = 4 * pi * self.dbt ** 2 / 4
                    Asy = 4 * pi * self.dbt ** 2 / 4
                else:
                    raise ValueError(f"Unknown lat_config ({self.lat_config})")
                 
                # flx and fly = lateral confining stress in x and y directions
                flx = ke * Asx / (self.s * dc) * self.fyt
                fly = ke * Asy / (self.s * bc) * self.fyt
                
                # Confinement effect from Chang, G. A., and Mander, J. B. (1994). 
                # Seismic Energy Based Fatigue Damage Analysis of Bridge Columns: Part I - 
                # Evaluation of Seismic Capacity. National Center for Earthquake Engineering 
                # Research, Department of Civil Engineering, State University of New York at 
                # Buffalo, Buffalo, New York. 
                # Confinement Effect on Strength (Section 3.4.3)
                x_bar = (flx + fly) / (2 * self.fc)
                r = max(flx, fly) / min(flx, fly)
                A_parameter = 6.886 - (0.6069 + 17.275*r) * exp(-4.989*r)
                B_parameter = 4.5 / (5 / A_parameter * (0.9849 - 0.6306 * exp(-3.8939*r)) - 0.1) - 5
                k1 = A_parameter * (0.1 + 0.9 / (1 + B_parameter * x_bar))
                fcc = self.fc * (1 + k1 * x_bar)
                # Confinement Effect on Ductility (Section 3.4.4)
                k2 = 5 * k1
                eps_prime_cc = self.eps_c * (1 + k2 * x_bar)

                ops.uniaxialMaterial("Concrete04", cover_concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)
                ops.uniaxialMaterial("Concrete04", core_concrete_material_id, -fcc, -eps_prime_cc, -1.0, self.Ec)
                confinement = True

            elif conc_mat_type == "Concrete04_no_confinement":
                ops.uniaxialMaterial("Concrete04", concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)
                confinement = False

            elif conc_mat_type == "Concrete01_no_confinement":
                ops.uniaxialMaterial("Concrete01", concrete_material_id, -self.fc, -2*self.fc/self.Ec, 0, 0.006)
                confinement = False

            elif conc_mat_type == "ENT":
                ops.uniaxialMaterial('ENT', concrete_material_id, self.Ec)
                confinement = False

            elif conc_mat_type == "Elastic":
                ops.uniaxialMaterial('Elastic', concrete_material_id, self.Ec)
                confinement = False

            else:
                raise ValueError(f"Concrete material {conc_mat_type} not supported")
            # endregion

            # region Define fibers and patches
            if (axis is None) or (axis == '3d'): 
                for i in self.reinforcement:
                    for index, value in enumerate(i.coordinates[0]):
                        ops.fiber(i.coordinates[1][index], value, i.Ab, steel_material_id)
                        if confinement:
                            negative_area_material_id = core_concrete_material_id
                        else:
                            negative_area_material_id = concrete_material_id
                        ops.fiber(i.coordinates[1][index], value, -i.Ab, negative_area_material_id)

                H = self.conc_cross_section.H
                B = self.conc_cross_section.B
                cdb = (H - self.reinforcement[0].By) / 2 - self.reinforcement[0].db / 2
                if self.dbt is not None:
                    cdb = cdb - self.dbt / 2
                
                if confinement:
                    nfy_cover = ceil(cdb * nfy / H)
                    nfx_cover = ceil(cdb * nfx / B)
                    nfy_core = ceil((H - 2 * cdb) * nfy / H)
                    nfx_core = ceil((B - 2 * cdb) * nfx / B)

                    ops.patch('rect', cover_concrete_material_id, nfy_cover, nfx,
                              -H / 2, -B / 2, -H / 2 + cdb, B / 2)

                    ops.patch('rect', cover_concrete_material_id, nfy_core, nfx_cover,
                              -H / 2 + cdb, -B / 2, H / 2 - cdb, -B / 2 + cdb)

                    ops.patch('rect', cover_concrete_material_id, nfy_cover, nfx,
                              H / 2 - cdb, -B / 2, H / 2, B / 2)

                    ops.patch('rect', cover_concrete_material_id, nfy_core, nfx_cover,
                              -H / 2 + cdb, B / 2 - cdb, H / 2 - cdb, B / 2)

                    ops.patch('rect', core_concrete_material_id, nfy_core, nfx_core,
                              -H / 2 + cdb, -B / 2 + cdb, H / 2 - cdb, B / 2 - cdb)
                else:
                    ops.patch('rect', concrete_material_id, nfy, nfx, -H / 2, -B / 2, H / 2, B / 2)

            elif axis in ['x', 'y']:
                for i in self.reinforcement:
                    if axis == 'x':
                        rebar_coords = i.coordinates[1]
                        H = self.conc_cross_section.H
                        B = self.conc_cross_section.B
                        nfd = nfy
                    else:
                        rebar_coords = i.coordinates[0]
                        H = self.conc_cross_section.B
                        B = self.conc_cross_section.H
                        nfd = nfx

                    for index, value in enumerate(rebar_coords):
                        ops.fiber(float(value), 0, i.Ab, steel_material_id)
                        if confinement:
                            negative_area_material_id = core_concrete_material_id
                        else:
                            negative_area_material_id = concrete_material_id
                        ops.fiber(float(value), 0, -i.Ab, negative_area_material_id)

                cdb = (H - self.reinforcement[0].By) / 2 - self.reinforcement[0].db / 2
                if self.dbt is not None:
                    cdb = cdb - self.dbt / 2
                if confinement:
                    nfd_cover = ceil(cdb * nfd / H)
                    nfd_core = ceil((H - 2 * cdb) * nfd / H)

                    core_fiber_height = (H - 2 * cdb) / nfd_core
                    cover_fiber_height = cdb / nfd_cover

                    ops.layer('straight', core_concrete_material_id, nfd_core, core_fiber_height * B,
                              -H / 2 + cdb, 0, H / 2 - cdb, 0)
                    ops.layer('straight', cover_concrete_material_id, nfd_cover, cover_fiber_height * B,
                              -H / 2, 0, -H / 2 + cdb, 0)
                    ops.layer('straight', cover_concrete_material_id, nfd_cover, cover_fiber_height * B,
                              H / 2 - cdb, 0, H / 2, 0)
                else:
                    fiber_height = H / nfd
                    ops.layer('straight', concrete_material_id, nfd, fiber_height * B, -H / 2, 0, H / 2, 0)
        
        elif type(self.conc_cross_section).__name__ == 'Circle':
        
            # region Define Concrete Material
            if conc_mat_type == "Concrete04":
                # Defined based on Mander, J. B., Priestley, M. J. N., and Park, R. (1988).
                # “Theoretical Stress-Strain Model for Confined Concrete.” Journal of Structural
                # Engineering, ASCE, 114(8), 1804―1826.
                
                if type(self.reinforcement[0]).__name__ != 'ReinfCirc':
                    raise ValueError(f"Reinforcement type {type(self.reinforcement[0]).__name__} not supported for this section type")
                if self.reinforcement[0].xc != 0 or self.reinforcement[0].yc != 0:
                    raise ValueError(f"Reinforcing pattern must be centered")
                if self.dbt is None:
                    raise ValueError("dbt must be defined")
                if self.s is None:
                    raise ValueError("s must be defined")

                # ds = core diameter based on center line of perimeter hoop
                ds = 2*self.reinforcement[0].rc + self.reinforcement[0].db + self.dbt
                
                # Ac = area of core of section enclosed by the center line of the permiter hoops
                Ac = pi/4*ds*ds
                
                # ρcc = ratio of area of longitudinal reinforcement to area of core of section
                ρcc = self.reinforcement[0].Ab * self.reinforcement[0].num_bars / Ac
                
                # sp = clear vertical spacing between spiral or hoop bars
                sp = self.s - self.dbt
                
                # ke = confinement effectiveness coefficient (Equation 15)
                ke = (1 - sp/(2*ds))/(1-ρcc)
                
                # flx and fly = lateral confining stress in x and y directions
                ρs = 4*self.Abt/(ds*self.s)  # Equation 17
                fl = 0.5*ke*ρs*self.fyt  # Equation 19
                
                # fcc = confined concrete strength
                fcc = self.fc*(-1.254 + 2.254*sqrt(1+7.94*fl/self.fc) - 2*fl/self.fc)
                
                # Confinement Effect on Ductility (Section 3.4.4 of Chang and Mander 1994)
                k1 = (fcc-self.fc)/fl
                k2 = 5 * k1
                x_bar = (fl + fl) / (2 * self.fc)
                eps_prime_cc = self.eps_c * (1 + k2 * x_bar)

                ops.uniaxialMaterial("Concrete04", cover_concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)
                ops.uniaxialMaterial("Concrete04", core_concrete_material_id, -fcc, -eps_prime_cc, -1.0, self.Ec)
                confinement = True
    
            elif conc_mat_type == "Concrete04_no_confinement":
                ops.uniaxialMaterial("Concrete04", concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)
                confinement = False

            elif conc_mat_type == "Concrete01_no_confinement":
                ops.uniaxialMaterial("Concrete01", concrete_material_id, -self.fc, -2*self.fc/self.Ec, 0, 0.006)
                confinement = False

            elif conc_mat_type == "ENT":
                ops.uniaxialMaterial('ENT', 2, self.Ec)
                confinement = False
    
            elif conc_mat_type == "Elastic":
                ops.uniaxialMaterial('Elastic', 2, self.Ec)
                confinement = False
    
            else:
                raise ValueError(f"Concrete material {conc_mat_type} not supported")
            # endregion
    
            # region Define fibers and patches
            if (axis is None) or (axis == '3d'): 
                for i in self.reinforcement:
                    for index, value in enumerate(i.coordinates[0]):
                        ops.fiber(i.coordinates[1][index], value, i.Ab, steel_material_id)
                        if confinement:
                            negative_area_material_id = core_concrete_material_id
                        else:
                            negative_area_material_id = concrete_material_id
                        ops.fiber(i.coordinates[1][index], value, -i.Ab, negative_area_material_id)
        
                d  = self.conc_cross_section.diameter
                max_fiber_size = d/max(nfx,nfy)
                if confinement:
                    ds = 2*self.reinforcement[0].rc + self.reinforcement[0].db + self.dbt
                    # Core Concrete
                    nfr = ceil(0.125*ds/max_fiber_size)              
                    nfc = ceil(0.25*ds*pi/max_fiber_size)
                    ops.patch('circ', core_concrete_material_id, nfc, nfr, 0, 0, 0,        0.125*ds, 0, 360)
                    nfc = ceil(0.5*ds*pi/max_fiber_size)
                    ops.patch('circ', core_concrete_material_id, nfc, nfr, 0, 0, 0.125*ds, 0.250*ds, 0, 360)
                    nfc = ceil(0.75*ds*pi/max_fiber_size)
                    ops.patch('circ', core_concrete_material_id, nfc, nfr, 0, 0, 0.250*ds, 0.375*ds, 0, 360)
                    nfc = ceil(ds*pi/max_fiber_size)
                    ops.patch('circ', core_concrete_material_id, nfc, nfr, 0, 0, 0.375*ds, 0.500*ds, 0, 360)
                    # Cover Concrete
                    nfr = ceil(0.5*(d-ds)/max_fiber_size)
                    nfc = ceil(d*pi/max_fiber_size)
                    ops.patch('circ', cover_concrete_material_id, nfc, nfr, 0, 0, 0.5*ds, 0.5*d, 0, 360)
                else:
                    nfr = ceil(0.125*d/max_fiber_size)              
                    nfc = ceil(0.25*d*pi/max_fiber_size)
                    ops.patch('circ', concrete_material_id, nfc, nfr, 0, 0, 0,       0.125*d, 0, 360)
                    nfc = ceil(0.5*d*pi/max_fiber_size)
                    ops.patch('circ', concrete_material_id, nfc, nfr, 0, 0, 0.125*d, 0.250*d, 0, 360)
                    nfc = ceil(0.75*d*pi/max_fiber_size)
                    ops.patch('circ', concrete_material_id, nfc, nfr, 0, 0, 0.250*d, 0.375*d, 0, 360)
                    nfc = ceil(d*pi/max_fiber_size)
                    ops.patch('circ', concrete_material_id, nfc, nfr, 0, 0, 0.375*d, 0.500*d, 0, 360)
            elif axis in ['x', 'y']:
                for i in self.reinforcement:
                    if axis == 'x':
                        rebar_coords = i.coordinates[1]
                    else:
                        rebar_coords = i.coordinates[0]

                    for index, value in enumerate(i.coordinates[1]):
                        ops.fiber(value, 0, i.Ab, steel_material_id)
                        if confinement:
                            negative_area_material_id = core_concrete_material_id
                        else:
                            negative_area_material_id = concrete_material_id
                        ops.fiber(value, 0, -i.Ab, negative_area_material_id)            
        
                d  = self.conc_cross_section.diameter
                max_fiber_size = d/nfy
                if confinement:
                    ds = 2*self.reinforcement[0].rc + self.reinforcement[0].db + self.dbt
                    nf = ceil(0.5*ds/max_fiber_size)
                    circ_patch_2d(core_concrete_material_id, nf, ds)
                    nf = ceil(0.5*d/max_fiber_size)
                    circ_patch_2d(cover_concrete_material_id, nf, d, Di=ds)
                else:
                    nf = ceil(0.5*d/max_fiber_size)
                    circ_patch_2d( concrete_material_id, nf, d)
            else:
                raise ValueError(f'Unknown option for axis: {axis}')
            # endregion

        elif type(self.conc_cross_section).__name__ == 'Obround':

            # region Define Concrete Material
            if conc_mat_type == "Concrete04":
                # Defined based on Mander, J. B., Priestley, M. J. N., and Park, R. (1988).
                # “Theoretical Stress-Strain Model for Confined Concrete.” Journal of Structural
                # Engineering, ASCE, 114(8), 1804―1826.

                if type(self.reinforcement[0]).__name__ != 'ReinfIntersectingLoops':
                    raise ValueError(f"Reinforcement type {type(self.reinforcement[0]).__name__} not supported for this section type")
                if self.reinforcement[0].xc != 0 or self.reinforcement[0].yc != 0:
                    raise ValueError(f"Reinforcing pattern must be centered")
                if self.dbt is None:
                    raise ValueError("dbt must be defined")
                if self.s is None:
                    raise ValueError("s must be defined")

                # ds = core diameter based on center line of perimeter hoop
                ds = self.reinforcement[0].D + self.reinforcement[0].db + self.dbt

                # Ac = area of core of section enclosed by the center line of the perimeter hoops
                Ac = pi / 4 * ds * ds

                # ρcc = ratio of area of longitudinal reinforcement to area of core of section
                ρcc = self.reinforcement[0].Ab * self.reinforcement[0].num_bars/2 / Ac

                # sp = clear vertical spacing between spiral or hoop bars
                sp = self.s - self.dbt

                # ke = confinement effectiveness coefficient (Equation 15)
                ke = (1 - sp / (2 * ds)) / (1 - ρcc)

                # flx and fly = lateral confining stress in x and y directions
                ρs = 4 * self.Abt / ds * self.s  # Equation 17
                fl = 0.5 * ke * ρs * self.fyt  # Equation 19

                # fcc = confined concrete strength
                fcc = self.fc * (-1.254 + 2.254 * sqrt(1 + 7.94 * fl / self.fc) - 2 * fl / self.fc)

                # Confinement Effect on Ductility (Section 3.4.4 of Chang and Mander 1994)
                k1 = (fcc - self.fc) / fl
                k2 = 5 * k1
                x_bar = (fl + fl) / (2 * self.fc)
                eps_prime_cc = self.eps_c * (1 + k2 * x_bar)

                ops.uniaxialMaterial("Concrete04", cover_concrete_material_id, -self.fc, -self.eps_c,
                                     -2 * self.eps_c, self.Ec)
                ops.uniaxialMaterial("Concrete04", core_concrete_material_id, -fcc, -eps_prime_cc,
                                     - 2 * eps_prime_cc, self.Ec)
                confinement = True
                
            elif conc_mat_type == "Concrete04_no_confinement":
                ops.uniaxialMaterial("Concrete04", 2, -self.fc, -self.eps_c, -2 * self.eps_c, self.Ec)
                confinement = False

            elif conc_mat_type == "Concrete01_no_confinement":
                ops.uniaxialMaterial("Concrete01", concrete_material_id, -self.fc, -2 * self.fc / self.Ec, 0, 0.006)
                confinement = False

            elif conc_mat_type == "ENT":
                ops.uniaxialMaterial('ENT', 2, self.Ec)
                confinement = False

            elif conc_mat_type == "Elastic":
                ops.uniaxialMaterial('Elastic', 2, self.Ec)
                confinement = False

            else:
                raise ValueError(f"Concrete material {conc_mat_type} not supported")
            # endregion

            # region Define fibers and patches
            if (axis is None) or (axis == '3d'):
                raise ValueError(f'3D option not supported for obround cross-sections yet')

            elif axis == 'x' or axis == 'y':

                if confinement:
                    negative_area_material_id = core_concrete_material_id
                else:
                    negative_area_material_id = concrete_material_id

                for i in self.reinforcement:
                    if axis == 'x':
                        for index, value in enumerate(i.coordinates[1]):
                            ops.fiber(value, 0, i.Ab, steel_material_id)
                            ops.fiber(value, 0, -i.Ab, negative_area_material_id)
                    elif axis == 'y':
                        for index, value in enumerate(i.coordinates[0]):
                            ops.fiber(value, 0, i.Ab, steel_material_id)
                            ops.fiber(value, 0, -i.Ab, negative_area_material_id)

                if axis == 'x':
                    nf = nfx
                elif axis == 'y':
                    nf = nfy

                if confinement:
                    obround_patch_2d_confined(cover_concrete_material_id, core_concrete_material_id, nf,
                                              self.conc_cross_section.D, self.conc_cross_section.a, ds,
                                              axis=axis)
                else:
                    obround_patch_2d(concrete_material_id, nf, 
                                     self.conc_cross_section.D, self.conc_cross_section.a,
                                     axis=axis)

            else:
                raise ValueError(f'Unknown option for axis: {axis}')
            # endregion

        else:
            raise ValueError(f"Concrete cross section {self.conc_cross_section.section_type} not supported")


def run_example():
    from libdenavit.section import Rectangle, ReinfRect

    # Define concrete cross section
    H = 40
    B = 20
    conc_cross_section = Rectangle(H, B)

    # Define reinforcement
    rhosr = 0.06
    nbB = 2
    nbH = 3
    cover = 0.15 * H
    Ab = H * B * rhosr / (2 * nbB + 2 * nbH - 4)
    reinforcement = ReinfRect(B - 2 * cover, H - 2 * cover, nbB, nbH, Ab)

    # Define materials
    fc = 4
    fy = 60
    units = 'US'

    # Define RC object
    section = RC(conc_cross_section, reinforcement, fc, fy, units)
    
    # Plot Section
    section.plot_section(show=False)

    # Plot Interaction Diagram
    angle = 0
    num_points = 40
    P, M, et = section.section_interaction_2d("x",num_points)
    ϕ = section.phi(et)
   
    plt.figure()
    plt.plot(M, -P, '-o', label="$M_x$ (nominal)")
    plt.plot(ϕ*M, -ϕ*P, '-s', label="$M_x$ (design)")
    plt.xlabel('Bending Moment (kip-in.)')
    plt.ylabel('Axial Compression (kips)')
    plt.legend()
    plt.show()


def run_example_2():
    D = 36
    a = 18
    Dc = D - 6
    nb = 20
    Ab = 1

    from libdenavit.section import Obround, ReinfIntersectingLoops
    conc_cross_section = Obround(D, a)

    # Define reinforcement
    reinforcement = ReinfIntersectingLoops(Dc, a, nb, Ab)

    # Define materials
    fc = 4
    fy = 60
    units = 'US'
    axis = 'y'
    # Define RC object
    section = RC(conc_cross_section, reinforcement, fc, fy, units)

    # Plot Section
    section.plot_section(show=False)

    # Plot Interaction Diagram
    angle = 0
    num_points = 40
    P, M, et = section.section_interaction_2d('x', num_points)
    ϕ = section.phi(et)

    plt.figure()
    plt.plot(M, -P, '-o', label="$M_x$ (nominal)")
    plt.plot(ϕ * M, -ϕ * P, '-s', label="$M_x$ (design)")
    plt.xlabel('Bending Moment (kip-in.)')
    plt.ylabel('Axial Compression (kips)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    # run_example()
    run_example_2()
