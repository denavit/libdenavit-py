import dataclasses
from math import pi
from libdenavit.design import available_strength

d = {
    '1/2':   0.5,
    '5/8':   0.625,
    '3/4':   0.75,
    '7/8':   0.875,
    '1':     1.0,
    '1-1/8': 1.125,
    '1-1/4': 1.25,
    '1-3/8': 1.375,
    '1-1/2': 1.5
}

dh_STD = {
    '1/2':   9/16,
    '5/8':   11/16,
    '3/4':   13/16,
    '7/8':   15/16,
    '1':     1.125,
    '1-1/8': 1.25,
    '1-1/4': 1.375,
    '1-3/8': 1.5,
    '1-1/2': 1.625
}

dh_OVS = {
    '1/2':   5/8,
    '5/8':   13/16,
    '3/4':   15/16,
    '7/8':   1.0625,
    '1':     1.25,
    '1-1/8': 1.4375,
    '1-1/4': 1.5625,
    '1-3/8': 1.6875,
    '1-1/2': 1.8125
}

Tb_GroupA = {
    '1/2':    12,
    '5/8':    19,
    '3/4':    28,
    '7/8':    39,
    '1':      51,
    '1-1/8':  64,
    '1-1/4':  81,
    '1-3/8':  97,
    '1-1/2': 118
}

Tb_GroupB = {
    '1/2':    15,
    '5/8':    24,
    '3/4':    35,
    '7/8':    49,
    '1':      64,
    '1-1/8':  80,
    '1-1/4': 102,
    '1-3/8': 121,
    '1-1/2': 148
}

Tb_GroupC = {
    '1':      90,
    '1-1/8': 113,
    '1-1/4': 143
}

Fnt = {
    'A307':      45,
    'GroupA-N':  90,
    'GroupA-X':  90,
    'GroupB-N': 113,
    'GroupB-X': 113,
    'GroupC-N': 150,
    'GroupC-X': 150
}

Fnv = {
    'A307':      27,
    'GroupA-N':  54,
    'GroupA-X':  68,
    'GroupB-N':  68,
    'GroupB-X':  84,
    'GroupC-N':  90,
    'GroupC-X': 113
}

mu = {
    'ClassA': 0.30,
    'ClassB': 0.50
}

@dataclasses.dataclass
class Bolt:
    """
    Parameters
    ----------
    d_str : str
        fractional bolt diameter
    bolt_type : str
        type of bolt
    hole_type : str, optional
        type of hole, default = STD
    slip_critical: bool
        is the connection slip-critical, default = False
    surface_type : str, optional
        type of surface
    """
    d_str: str
    bolt_type: str
    hole_type: str = 'STD'
    slip_critical: bool = False
    surface_type: str = None
    
    strength_type = 'design'
    deformation_considered = True  # deformation at the bolt hole at service load is a design consideration
    Du = 1.13
    hf = 1.0
        
    @property
    def d(self):
        return d[self.d_str]

    @property
    def Ab(self):
        Ab = 0.25*pi*self.d**2
        return Ab
        
    @property
    def dh(self):
        if self.hole_type == 'STD':
            dh = dh_STD[self.d_str]
        elif self.hole_type == 'OVS':
            dh = dh_OVS[self.d_str]
        else:
            raise Exception('Unknown hole_type: %s' % self.hole_type)
        return dh
     
    @property
    def Fnt(self):
        return Fnt[self.bolt_type] 
     
    @property
    def Fnv(self):
        return Fnv[self.bolt_type]       

    @property
    def Tb(self):
        if self.bolt_type == 'GroupA-N' or self.bolt_type == 'GroupA-X':
            return Tb_GroupA[self.d_str]
        if self.bolt_type == 'GroupB-N' or self.bolt_type == 'GroupB-X':
            return Tb_GroupB[self.d_str]
        if self.bolt_type == 'GroupC-N' or self.bolt_type == 'GroupC-X':
            return Tb_GroupC[self.d_str]
        else:
            raise Exception('Unknown bolt_type: %s' % self.bolt_type)  

    @property
    def mu(self):
        return mu[self.surface_type] 

    def rn_bolt_tension(self,frv=0.0):
        if frv == 0.0:
            rn = self.Fnt*self.Ab
        else:
            Fnv = available_strength(self.Fnv,self.strength_type,0.75,2.00)
            rn = min(1.3*self.Fnt-frv*self.Fnt/Fnv,self.Fnt)*self.Ab
        return available_strength(rn,self.strength_type,0.75,2.00)

    def rn_bolt_shear(self,ns):
        rn = self.Fnv*self.Ab*ns
        return available_strength(rn,self.strength_type,0.75,2.00)
        
    def rn_bearing(self,t,Fu):
        if self.deformation_considered:
            rn = 2.4*self.d*t*Fu
        else:
            rn = 3.0*self.d*t*Fu
        return available_strength(rn,self.strength_type,0.75,2.00)

    def rn_tearout(self,lc,t,Fu):
        if self.deformation_considered:
            rn = 1.2*lc*t*Fu
        else:
            rn = 1.5*lc*t*Fu
        return available_strength(rn,self.strength_type,0.75,2.00)   
        
    def rn_slip(self,ns):
        rn = self.mu*self.Du*self.hf*self.Tb*ns        
        if self.hole_type == 'STD':
            return available_strength(rn,self.strength_type,1.00,1.50)
        elif self.hole_type == 'OVS':
            return available_strength(rn,self.strength_type,0.85,1.76)
        else:
            raise Exception('Unknown hole_type: %s' % self.hole_type)
    
def run_example():    
    b = Bolt('3/4','GroupA-N',hole_type='OVS',surface_type='ClassA')

    print(f'd = {b.d} in.')
    print(f'dh = {b.dh} in.')
    print(f'Ab = {b.Ab:.3f} in.^2')
    print(f'φrn = {b.rn_bolt_shear(2):.3f} kips (shear rupture, double shear)')
    print(f'φrn = {b.rn_slip(2):.3f} kips (slip, 2 slip planes)')
    
if __name__ == "__main__":
    run_example()