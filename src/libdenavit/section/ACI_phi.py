import numpy as np

def ACI_phi(transverse_reinf_type, et, ety):
    
    # Table 21.2.2, ACI 318-19 
    if transverse_reinf_type.lower() in ['spiral', 'spiral']:
        phicc = 0.75
    elif transverse_reinf_type.lower() in ['other', 'ties']:
        phicc = 0.65
    else:
        raise ValueError(f'unknown type {transverse_reinf_type}')

    if type(ety) == str and ety.lower() == 'grade60':
        ety = 0.002

    phitc = 0.9

    def compute_phi(et):
        # Tension Controlled
        if et >= ety + 0.003:
            phi = phitc

        # Compression Controlled
        elif et <= ety:
            phi = phicc

        # Transition
        else:
            phi = phicc + (phitc - phicc) * (et - ety) / 0.003
        
        return phi

    if isinstance(et,np.ndarray):
        phi = np.empty_like(et)
        for i,iet in enumerate(et):
            phi[i] = compute_phi(iet)
        return phi
        
    phi = compute_phi(et)   

    return phi
