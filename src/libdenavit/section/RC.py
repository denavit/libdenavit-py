from math import sqrt, pi, ceil
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberSection, ACI_phi
import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops


class RC:
    transverse_reinf_type = ''
    treat_reinforcement_as_point = True
    _reinforcement = None
    _Ec = None
    _Es = None
    _eps_c = None
    _Abt = None

    def __init__(self, conc_cross_section, reinforcement, fc, fy, units, dbt=None, s=None, fyt=None, lat_config=None):
        self.conc_cross_section = conc_cross_section
        self.reinforcement = reinforcement
        self.fc = fc
        self.fy = fy
        self.units = units
        self.dbt = dbt
        self.s = s
        self.fyt = fyt
        self.lat_config = lat_config

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
        if self._Abt is None:
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
        p0 = 0.85 * self.fc * (self.Ag - self.Asr) + self.fy * self.Asr
        return p0

    @property
    def pnco(self):
        # See Section 22.4.2 of ACI318-19
        if self.transverse_reinf_type.lower() == 'ties':
            pnco = 0.80 * self.p0
        elif self.transverse_reinf_type.lower() == 'spiral' or 'spirals':
            pnco = 0.85 * self.p0
        else:
            raise ValueError("Unknown transverse_reinf_type")
        return pnco

    def phi(self, et):
        f = ACI_phi(self.transverse_reinf_type, et, self.fy / self.Es)
        return f

    def plot_section(self, line_width=2, save_name=None):
        plt.figure()
        self.conc_cross_section.plot_section(line_width)
        for i in range(len(self.reinforcement)):
            self.reinforcement[i].plot_section()
        plt.box(False)
        plt.axis('scaled')

        if save_name:
            plt.savefig(save_name)
        else:
            plt.show()

        plt.clf()
        plt.close()

    def aci_strain_compatibility_object(self):
        id_conc = 1
        id_reinf = 2
        fs = self.fiber_section_object(id_conc, id_reinf)
        fs.print_section_properties()
        scACI = AciStrainCompatibility(fs)
        # add concrete boundaries
        x, y, r = self.conc_cross_section.boundary_points
        for i in range(len(x)):
            scACI.add_concrete_boundary(x[i], y[i], r[i])
        # add steel boundary
        for i in self.reinforcement:
            x, y = i.coordinates

            for j in range(len(x)):
                if self.treat_reinforcement_as_point:
                    scACI.add_steel_boundary(x[j], y[j], 0)
                else:
                    raise ValueError("Not implemented yet")

        scACI.maxCompressiveStrength = -self.pnco
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.fy, self.Es)
        return scACI

    def fiber_section_object(self, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        self.conc_cross_section.add_to_fiber_section(fs, id_conc)
        for i in self.reinforcement:
            i.add_to_fiber_section(fs, id_reinf, id_conc)
        return fs

    def section_interaction_2d(self):
        scACI = self.aci_strain_compatibility_object()

        plot_P = []
        plot_Mx = []
        plot_et = []
        scACI.build_data()

        upper = self.conc_cross_section.H * 1.51
        lower = -self.conc_cross_section.H * 2
        r1 = np.array([-1e99])
        r2 = np.arange(lower, lower / 2, 0.2)
        r3 = np.arange(lower / 2, upper, 0.2)
        r4 = np.array([1e99])
        d = np.concatenate((r1, r2, r3, r4))

        for i in d:
            P, Mx, My, et = scACI.compute_point(0, i, 0)
            plot_P.append(-P)
            plot_Mx.append(Mx)
            plot_et.append(et)

        return plot_P, plot_Mx, plot_et

    def plot_interaction_diagram(self, save_name=None):
        plot_P, plot_Mx, plot_et = self.section_interaction_2d()
        phi = []
        for i in plot_et:
            phi.append(self.phi(i))
        plt.plot(plot_Mx, plot_P)
        phi = np.array(phi)
        plot_Mx = np.array(plot_Mx)
        plot_P = np.array(plot_P)
        plt.plot(phi * plot_Mx, phi * plot_P)
        plt.xlabel('M')
        plt.ylabel('P')
        plt.title("Interaction Diagram")

        if save_name:
            plt.savefig(f'{save_name}')
        else:
            plt.show()

        plt.clf()
        plt.close()

    def build_ops_fiber_section(self, section_id, start_material_id, steel_mat_type, conc_mat_type, nfy, nfx, GJ=1.0e6):
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
            ops.uniaxialMaterial("ReinforcingSteel", steel_material_id, self.fy, self.fy * 1.5, self.Es, self.Es / 2, 0.002, 0.008)

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
                    raise ValueError(f"Reinforcing patter must be centered")
                
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
                flx = Asx / (self.s * dc) * self.fyt
                fly = Asy / (self.s * bc) * self.fyt
                
                # Chang & Mander's 1994 concrete model- section 3.4.3-4
                q = max(flx, fly) / min(flx, fly)
                x_prime = (flx + fly) / (2 * self.fc)
                A_parameter = 6.886 - (0.6069 + 17.275*q) * np.e ** (-4.989*q)
                B_parameter = 4.5 / (5 / A_parameter * (0.9849 - 0.6306* np.e ** (-3.8939*q)) - 0.1) - 5
                k1 = A_parameter * (0.1 + 0.9 / (1 + B_parameter * x_prime))
                fcc = self.fc * (1 + k1 * x_prime)

                k2 = 5 * k1
                eps_prime_cc = self.eps_c * (1 + k2 * x_prime)

                ops.uniaxialMaterial("Concrete04", cover_concrete_material_id, -self.fc, -self.eps_c, -2 * self.eps_c, self.Ec)
                ops.uniaxialMaterial("Concrete04", core_concrete_material_id, -fcc, -eps_prime_cc, - 2 * eps_prime_cc, self.Ec)
                confinement = True

            elif conc_mat_type == "Concrete04_no_confinement":
                ops.uniaxialMaterial("Concrete04", concrete_material_id, -self.fc, -self.eps_c, -2 * self.eps_c, self.Ec)
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
            for i in self.reinforcement:
                for index, value in enumerate(i.coordinates[0]):
                    ops.fiber(value, i.coordinates[1][index], i.Ab, steel_material_id)
                    if confinement:
                        negative_area_material_id = core_concrete_material_id
                    else:
                        negative_area_material_id = concrete_material_id
                    ops.fiber(value, i.coordinates[1][index], -i.Ab, negative_area_material_id)

            H = self.conc_cross_section.H
            B = self.conc_cross_section.B
            cdb = (H - self.reinforcement[0].Bx)/2 - self.reinforcement[0].db / 2
            if self.dbt is not None:
                cdb = cdb - self.dbt / 2
            
            if confinement:
                ops.patch('rect', cover_concrete_material_id, ceil(cdb * nfy / H), nfx, 
                    -H/2, -B/2, -H/2 + cdb, B/2)
                
                ops.patch('rect', cover_concrete_material_id, ceil((H - 2 * cdb) * nfy / H), ceil(cdb * nfx / B), 
                    -H/2 + cdb, -B/2, H/2 - cdb, -B/2 + cdb)
                
                ops.patch('rect', cover_concrete_material_id, ceil(cdb * nfy / H), nfx, 
                    H/2 - cdb, -B/2, H/2, B/2)
                
                ops.patch('rect', cover_concrete_material_id, ceil((H - 2 * cdb) * nfy / H), ceil(cdb * nfx / B),
                    -H/2 + cdb, B/2 - cdb, H/2 - cdb, B/2)
                
                ops.patch('rect', core_concrete_material_id, ceil((H - 2 * cdb) * nfy / H), ceil((B - 2 * cdb) * nfx / B),
                    -H/2 + cdb, -B/2 + cdb, H/2 - cdb, B/2 - cdb)
            else:
                ops.patch('rect', concrete_material_id, nfy, nfx, -H/2, -B/2, H/2, B/2)
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
    nbH = 2
    cover = 0.15 * H
    Ab = H * B * rhosr / (2 * nbB + 2 * nbH - 4)
    reinforcement = ReinfRect(B - 2 * cover, H - 2 * cover, nbB, nbH, Ab)

    # Define materials
    fc = 4
    fy = 60
    units = 'US'

    # Define RC object
    section = RC(conc_cross_section, reinforcement, fc, fy, units)
    section.transverse_reinf_type = 'ties'
    section.plot_section()
    section.plot_interaction_diagram()


if __name__ == '__main__':
    run_example()
