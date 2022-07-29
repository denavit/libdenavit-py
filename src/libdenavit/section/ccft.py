from math import pi
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberCirclePatch, FiberSection
from libdenavit.section.database import reinforcing_bar_database
import numpy as np


class CCFT:
    _Es = None
    _Ec = None
    
    def __init__(self, D, t, Fy, fc, units, num_bars=0, bar_size=None, Ab=None, db=None, Fylr=None, Dp=None, neglect_local_buckling=False):
        # Main Parameters
        self.D = D
        self.t = t
        self.Fy = Fy
        self.fc = fc
        self.units = units
        
        # Reinforcing
        self.num_bars = num_bars
        if bar_size is None:
            self.Ab = Ab
            self.db = db
        else:
            bar_data = reinforcing_bar_database[bar_size]
            self.Ab = bar_data['area']
            self.db = bar_data['diameter']
        self.Fylr = Fylr    # Yield strength of longitudinal reinforcing
        self.Dp = Dp        # Distance from outside surface of steel tube to center of longitudinal reinforcing
        
        # Options
        self.neglect_local_buckling = neglect_local_buckling

    @property
    def Es(self):
        if self._Es is not None:
            return self._Es

        if self.units.lower() == 'us':
            return 29000

        if self.units.lower() == 'si':
            return 200000

        raise ValueError(f'Es is not set and default is not implemented for {self.units = }')

    @Es.setter
    def Es(self, x):
        self._Es = x
         
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
    def As(self):
        return 0.25*pi*(self.D**2 - (self.D-2*self.t)**2)
        
    @property
    def Asr(self):
        if self.num_bars == 0:
            return 0.0
        else:
            return self.num_bars*self.Ab
    
    @property
    def Ac(self):
        return 0.25*pi*(self.D-2*self.t)**2 - self.Asr
    
    def reinforcing_coordinates(self):
        angles = np.linspace(0,2*pi,self.num_bars,endpoint=False)
        r = 0.5*self.D - self.Dp
        x = r*np.sin(angles)
        y = r*np.cos(angles)
        return x,y
    
    def aci_strain_compatibility_object(self):
        id_steel = 1
        id_conc = 2
        id_reinf = 3       
        fs = self.fiber_section_object(id_steel, id_conc, id_reinf)
        scACI = AciStrainCompatibility(fs)
        scACI.add_concrete_boundary(0.0, 0.0, 0.5*self.D-self.t)
        scACI.add_steel_boundary(0.0, 0.0, 0.5*self.D)
        scACI.add_material(id_steel, 'steel', self.Fy, self.Es)
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.Fylr, self.Es)
        return scACI    
    
    def fiber_section_object(self, id_steel, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        # Steel
        fs.add_fibers(FiberCirclePatch(0,0,0.5*self.D-self.t,0.5*self.D,id_steel))
        # Concrete 
        fs.add_fibers(FiberCirclePatch(0,0,0,0.5*self.D-self.t,id_conc))
        # Reinforcement
        if self.num_bars > 0:
            x,y = self.reinforcing_coordinates()
            for i in range(self.num_bars):
                fs.add_fibers(FiberSingle(self.Ab,x[i],y[i],id_reinf,id_conc))
        return fs
 


def run_example():
    #section = CCFT(24,0.25,50,6,'US')
    section = CCFT(24,0.25,50,6,'US',num_bars=8, bar_size='#8', Fylr=60, Dp=3.0)
    
    fs = section.fiber_section_object(1,2,3)
    fs.print_section_properties()
    
    print(f'As  = {section.As}')
    print(f'Ac  = {section.Ac}')
    print(f'Asr = {section.Asr}')
    
    
if __name__ == '__main__':
    run_example()   
    
    