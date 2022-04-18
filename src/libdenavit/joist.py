from math import floor, ceil
import numpy as np
import libdenavit

class BaseJoistClass:
    def moment_strength_ratio(self,x,Mr):
        (positive_envelope,negative_envelope) = self.moment_strength_envelope(x)
        
        strength_ratio = np.zeros_like(Mr)
        for i in range(len(Mr)):
            if abs(Mr[i]) <= 1e-8:
                strength_ratio[i] = 0.0
            elif Mr[i] < 0.0:
                strength_ratio[i] = Mr[i]/negative_envelope[i]
            else:
                strength_ratio[i] = Mr[i]/positive_envelope[i]
                
        return strength_ratio
        
        
    def shear_strength_ratio(self,x,Vr):
        (positive_envelope,negative_envelope) = self.shear_strength_envelope(x)
        
        strength_ratio = np.zeros_like(Vr)
        for i in range(len(Vr)):
            if Vr[i] < 0.0:
                strength_ratio[i] = Vr[i]/negative_envelope[i]
            else:
                strength_ratio[i] = Vr[i]/positive_envelope[i]
        
        return strength_ratio
        
    def max_strength_ratio(self,x,Mr,Vr):
        strength_ratio = max(
            max(self.moment_strength_ratio(x,Mr)),
            max(self.shear_strength_ratio(x,Vr)))
        return strength_ratio
    

class OpenWebSteelJoist(BaseJoistClass):
    
    def __init__(self,strength_type,span_ft,wTL_plf,wLL_plf):
        self.strength_type = strength_type
        self.span_ft = span_ft
        self.wTL_plf = wTL_plf # Total load that can be supported
        self.wLL_plf = wLL_plf # Load that causes L/360 deflection
        self.minimum_shear_reversal_strength_ratio = 0.0
        self.use_proposed_shear_strength = True
        
    def moment_strength_envelope(self,x,units='lb-ft'):
        
        # NOTE: calculations are performed as if the design length equals the span. 
        L = self.span_ft
        
        # Compute envelopes
        positive_envelope = 0.5*self.wTL_plf*x*(L-x)
        negative_envelope = np.zeros_like(x)
        
        # Convert to desired units
        ucf = libdenavit.unit_conversion_factor('lb-ft', units)               
        return (ucf*positive_envelope,ucf*negative_envelope)

    def shear_strength_envelope(self,x,units='lbf',spread_dist=1e-6):
        
        # Create copy of x and spread out pairs of values
        x2 = x.copy()
        spread_out_pairs(x2,spread_dist)
        
        # Compute envelopes
        positive_envelope = np.zeros_like(x2)
        negative_envelope = np.zeros_like(x2)
        
        # NOTE: calculations are performed as if the design length equals the span. 
        L = self.span_ft
        R = 0.5*self.wTL_plf*L # End reaction
        
        for i in range(len(x2)):
        
            # Basic shear strength
            if x2[i] <= 0.5*self.span_ft:
                positive_envelope[i] = max(0.25*R,self.wTL_plf*(0.5*L-x2[i]))
                negative_envelope[i] = -self.minimum_shear_reversal_strength_ratio*R
            else:
                positive_envelope[i] = self.minimum_shear_reversal_strength_ratio*R
                negative_envelope[i] = min(-0.25*R,self.wTL_plf*(0.5*L-x2[i]))           

            if self.use_proposed_shear_strength:
                # Use minimum shear strength provisions that are not in the 2020 SJI Specification, but
                # that are required of all SJI members
                positive_envelope[i] = max(positive_envelope[i],0.75*R*(1-x2[i]/L)**2 + 0.25*R*(1-2*x2[i]/L))
                negative_envelope[i] = min(negative_envelope[i],-0.75*R*(x2[i]/L)**2 - 0.25*R*(2*x2[i]/L-1))

        # Convert to desired units
        ucf = libdenavit.unit_conversion_factor('lbf', units)
        return (ucf*positive_envelope,ucf*negative_envelope)
    
    def moment_of_inertia(self):
        I = 26.767e-6*self.wLL_plf*(self.span_ft-1/3)**3
        return I
        
 
class JoistGirder(BaseJoistClass):
   
    def __init__(self,strength_type,span_ft,d_in,num_spaces,P_kips):
        self.strength_type = strength_type
        self.span_ft = span_ft
        self.d_in = d_in
        self.num_spaces = num_spaces
        self.P_kips = P_kips
        
    def moment_strength_envelope(self,x,units='lb-ft'):
        
        # NOTE: calculations are performed as if the design length equals the span. 
        L = self.span_ft
        
        # Compute envelopes
        positive_envelope = np.zeros_like(x)
        negative_envelope = np.zeros_like(x)
        for i in range(len(x)):
            s1 = floor(self.num_spaces*x[i]/L)
            M1 = self.P_kips*L*(self.num_spaces-s1)*s1/(2*self.num_spaces)
            M2 = self.P_kips*L*(self.num_spaces-(s1+1))*(s1+1)/(2*self.num_spaces)
            x1 = s1/self.num_spaces*L
            x2 = (s1+1)/self.num_spaces*L
            positive_envelope[i] = M1 + (M2-M1)*(x[i]-x1)/(x2-x1)
        
        # Convert to desired units
        ucf = libdenavit.unit_conversion_factor('kip-ft', units)             
        return (ucf*positive_envelope,ucf*negative_envelope)

    def shear_strength_envelope(self,x,units='lbf',spread_dist=1e-6):
        
        # Create copy of x and spread out pairs of values
        x2 = x.copy()
        spread_out_pairs(x2,spread_dist)
        
        # Compute envelopes
        positive_envelope = np.zeros_like(x2)
        negative_envelope = np.zeros_like(x2)
        
        # NOTE: calculations are performed as if the design length equals the span. 
        L = self.span_ft
        R = 0.5*self.P_kips*(self.num_spaces-1) # End reaction
        
        for i in range(len(x2)):
            if x2[i] == L:
                s1 = self.num_spaces
            else:
                s1 = floor(self.num_spaces*x2[i]/L) + 1
            Vi = self.P_kips*((self.num_spaces+1)/2 - s1) 
            
            if x2[i] <= 0.5*self.span_ft:
                positive_envelope[i] = max(Vi,0.25*R)
                negative_envelope[i] = -0.25*positive_envelope[i]               
            else:
                negative_envelope[i] = min(Vi,-0.25*R)
                positive_envelope[i] = -0.25*negative_envelope[i]   
                
        # Convert to desired units
        ucf = libdenavit.unit_conversion_factor('kips', units)
        return (ucf*positive_envelope,ucf*negative_envelope)
   
    def moment_of_inertia(self):
        if self.strength_type == 'ASD':
            I = 0.027*self.num_spaces*self.P_kips*(self.span_ft-1/3)*self.d_in
        elif self.strength_type == 'LRFD':
            I = 0.018*self.num_spaces*self.P_kips*(self.span_ft-1/3)*self.d_in
        else:
            raise ValueError(f'Unknown strength type: {self.strength_type}')
        return I


def spread_out_pairs(x,spread_dist):
    # Spread out pairs of points
    len_x = len(x)
    i = 0
    while i < (len_x-1):
        if x[i] == x[i+1]:
            # Pair found
            if i+2 < len_x:
                # Pair is in the middle of array
                if x[i] == x[i+2]:
                    # More than a pair
                    # don't adjust values and incrmeent past
                    while x[i] == x[i+1]:
                        i += 1
                else:
                    x[i]   -= spread_dist
                    x[i+1] += spread_dist
                    i += 2    
            else:
                # Pair is at end of array
                x[i]   -= spread_dist
                x[i+1] += spread_dist
                i += 2
        else:
            # No pair
            i += 1
    return
